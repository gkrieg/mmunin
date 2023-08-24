#!/usr/bin/python

import sys
from optparse import OptionParser
import pickle as pkl

from ilp import find_shortest_hyperpath_parallel 
from pathheuristic import get_doubly_reachable_graph, tail_path_heuristic
from cutfinder import convert_vertex_cuts_to_edge_cuts

from halp.directed_hypergraph import DirectedHypergraph

## GLOBAL VARIABLES
ROOTDIR = '/groups/kece/skrieger/mmunin'

def main(args):
    opts = parseOptions(args)
    outputfile = f'results/{opts.target[0][-8:]}results.txt'

    if opts.pathwayreconstruct:
        #read through hyperedges file picking out those corresponding to the ones from a pathway
        pathwayedges = []

        import requests
        pathwayrequest = requests.get('https://www.pathwaycommons.org/pc2/search.json?q=pathway:{}&datasource=reactome&type=conversion'.format(opts.pathway[0]))
        if 'json' in pathwayrequest.headers['Content-Type']:
            jsonobject = pathwayrequest.json()
            for edge in jsonobject['searchHit']:
                pathwayedges.append(edge['uri'])

        Hpathway,_,__ = make_hypergraph(ROOTDIR+'/hypergraphs/{0}'.format(opts.name),select_edges=pathwayedges)
        pathwayedges = Hpathway.get_hyperedge_id_set()
        sources = set([v for v in Hpathway.node_iterator() if len(Hpathway.get_backward_star(v)) == 0])
        targets = set([v for v in Hpathway.node_iterator() if len(Hpathway.get_forward_star(v)) == 0])
        node_dict = pkl.load(open('node_dict.pkl','rb'))
        print(len(sources),[node_dict[s] for s in sources])
        print(len(targets),[node_dict[t] for t in targets])
        print('number of hyperedges in the reactome pathway: ',len(pathwayedges))

    H,_,__ = make_hypergraph(ROOTDIR+'/hypergraphs/{0}'.format(opts.name))
    if not opts.pathwayreconstruct:
        print(('target: {0}'.format(opts.target)))
        a = H.get_backward_star(opts.target[0])
        H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name)
    elif opts.pathwayreconstruct:
        H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',sources,targets,keep_sources=True)
        source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,single_targets_edge=True)

    if opts.cyclicparallel or opts.pathwayreconstruct:
        ilpname = '%s/results/reactome-%s-hypergraph-cyclic.lp' % (ROOTDIR,opts.name)
        node_dict = pkl.load(open('node_dict.pkl','rb'))
        H2 = get_doubly_reachable_graph(H,source,target,node_dict)
        print(len(H2.get_hyperedge_id_set()))
        sourcecut = set([source])
        sinkcut = H2.get_node_set()
        if target not in H2.get_node_set():
            print('target not reachable from source')
        sinkcut.remove(target)
        vertexheadcuts = (sourcecut,sinkcut)
        edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H2,vertexheadcuts,node_dict)

        if opts.pathwayreconstruct:
            taildistancelist,Hnew,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
            reactomehit = []
            both = []
            ilponly = []
            for e in P:
                for f in pathwayedges:
                    if edgeequals(Hnew,e,Hpathway,f) == True:
                        both.append(e)
                        reactomehit.append(f)
                        break
                else:
                    ilponly.append(e)
            reactomeonly = []
            for e in pathwayedges:
                if e not in reactomehit:
                    reactomeonly.append(e)
            print('{} edges were in the reactome pathway originally'.format(len(pathwayedges)))
            print('{} edges were not recovered from the reactome pathway: {}'.format(len(reactomeonly),reactomeonly))
            print('{} edges were recovered from the reactome pathway: {}'.format(len(both),both))
            print('{} edges were outside the reactome pathway: {}'.format(len(ilponly),ilponly))
        else:
            find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,vertexheadcuts=vertexheadcuts,outputfile=outputfile)

def edgeequals(G,e,H,f):
    ehead = G.get_hyperedge_head(e)
    etail = G.get_hyperedge_tail(e)
    fhead = H.get_hyperedge_head(f)
    ftail = H.get_hyperedge_tail(f)
    if len(ehead) != len(fhead):
        return False
    if len(etail) != len(ftail):
        return False
    for v in ehead:
        if v not in fhead:
            return False
    for v in etail:
        if v not in ftail:
            return False
    return True


def make_hypergraph(file_prefix,delim=';',sep='\t',keep_singleton_nodes=False,target=None,select_edges=[]):
    hypernodes = {}
    with open(file_prefix+'-hypernodes.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split(sep)
            if len(row) == 1:
                hypernodes[row[0]] = ['OtherComplexes-FIX']
            else:
                hypernodes[row[0]] = row[1].split(delim)
    print(('%d hypernodes from hypernodes file' % (len(hypernodes))))
    identifier2id = {}
    id2identifier = {}
    H = DirectedHypergraph()
    if keep_singleton_nodes:
        for n in hypernodes:
            H.add_node(n)

    skipped1 = 0
    skipped2 = 0
    tailsizes = []
    headsizes = []
    selfloops = []
    noinselfloops = 0
    indegree = []
    outdegree = []
    numtargets = 0
    numsources = 0

    with open(file_prefix+'-hyperedges.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            #row = line.strip().split(sep)
            row = line.strip().split()
            tail = set()
            head = set()

            ## Tail includes tail and regulators.
            ## Head includes head.
            if row[0] != 'None' and row[0] != '':
                tail.update(row[0].split(delim))
            if row[1] != 'None' and row[1] != '':
                head.update(row[1].split(delim))
            if row[2] != 'None' and row[2] != '':
                tail.update(row[2].split(delim))
            #These are the negative regulators!
            #if row[3] != 'None':
                #tail.update(row[3].split(delim))
            hedge_id = row[4]

            ## THIS IS A HACK FOR NOW ( should be incorporated in the make-hypergraph.py code)
            ## IGnore any reactions that have a Reactome Identifier (e.g. has "HSA") instead of
            ## a PAthway Commons identifier.
            if any(['HSA' in s for s in tail]+['HSA' in s for s in head]):
                skipped1+=1
            elif len(tail)==0 or len(head)==0:
                skipped2+=1
            elif select_edges == [] or hedge_id in select_edges:
                hid = H.add_hyperedge(tail,head,identifier=hedge_id)
                tailsizes.append(len(tail))
                headsizes.append(len(head))
                intersection = tail.intersection(head)
                if len(intersection) > 0:
                    selfloops.append([v for v in intersection])

                identifier2id[hedge_id] = hid
                id2identifier[hid] = hedge_id

    print(('%d reactions skipped because of Reactome identifier' % (skipped1)))
    print(('%d reactions skipped because of an empty tail or head' % (skipped2)))
    ## annotate nodes
    num_hypernodes = 0
    for node in H.get_node_set():
        if node in hypernodes and hypernodes[node] != [node]:
            H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
            num_hypernodes+=1
        else:
            H.add_node(node,is_hypernode=False,hypernode_members=[])

        H.add_node(node)

    return H, identifier2id, id2identifier

        

def getSourcesTargets2(name,H,graph_type,source_list,target_list,keep_sources=False):
    if keep_sources == True:
        sources = set(source_list)
        sources.add('http://pathwaycommons.org/pc12/Complex_4d1995472d186b954709879cb288a668')
        sources.add('http://pathwaycommons.org/pc12/Protein_ba155426a839df372d91bd7b7413a4b4')
        sources.add('http://pathwaycommons.org/pc12/Complex_e00e57513597efe208923db55fc2da3f')
        sources.add('http://pathwaycommons.org/pc12/Protein_75fb1a098e0dd9795f76102b820ae5f6')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_d300fb612cb231cd91749edbb44281bf')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_274b04f4b237e66066509730fdd27a5e')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_1eef4b089daf104d1c8d5fbe4c04970e')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_2669fc5e07ba41847df68257aeb61939')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_e48664ebaf9f7fd98ae77fc2e8e80ebd')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_ea368c90d5a68cdfa5846928f1f714c1')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_738a4cedb1e48bbc1a7199a99ade5afa')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_98cda7817ff4c0f3daf6e1652ef1aa83')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_789b3871b86d95ed2880cdad43be982c')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_21d72a6d0423a7aa244a0316c29aec13')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_b4d4dc64c46d4c13974cedb37d4518c3')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_ea368c90d5a68cdfa5846928f1f714c1')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_738a4cedb1e48bbc1a7199a99ade5afa')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_98cda7817ff4c0f3daf6e1652ef1aa83')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_d2967b8379ca02a85352118fa7cca5f5')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_d314f1d40db097b480c9d002b450fd31')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_73353698b897a5cbcb452e7cbd37cdc5')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_21903e848f1318bc56ad8d75b5d70796')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_1a0b828dcd30fd0f85cbec480fe3dad7')
        sources.add('http://pathwaycommons.org/pc12/Protein_268dd140e5cc229f5c2accd006e19d85')
        sources.add('http://pathwaycommons.org/pc12/Complex_3241cfa7cf77f2454be4bbc745771fa5')
        sources.add('http://pathwaycommons.org/pc12/Complex_13a235ef8f2a7f846c53663cedb2e5e5')
        sources.add('http://pathwaycommons.org/pc12/Complex_cf2e225345ac7efc8499273bfdf6d23d')
        sources.add('http://pathwaycommons.org/pc12/Complex_5ab9537e0ff27689747f16b1772840fd')
        sources.add('http://pathwaycommons.org/pc12/Complex_13b738018ef2d0bbcc20d246ed183fd2')
        sources.add('http://pathwaycommons.org/pc12/Complex_0c19e31fff33d3d80c186b60abe1c1e3')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_3cf16238f3a29742331da93a6397eab8')
        sources.add('http://pathwaycommons.org/pc12/Complex_8635539c551e9f5c2e90b35ff40e221b')
        sources.add('http://pathwaycommons.org/pc12/SmallMolecule_5c7e3625bcca8ef4e7e1e09932cd0e5f')
        sources.add('http://pathwaycommons.org/pc12/Complex_a8d85da8ff183e6ea3ac1e74ff81aa87')
        #sources.add('')
        targets = set(target_list)
    if name == 'allpid':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Complex_81ba3b0707b6c6abd477dd59315147f4'])
        if ':' in target_list[0]:
            targets = set([target_list[0]])
    elif name == 'allreactome':
        sources = set(['http://pathwaycommons.org/pc12/Complex_eb845eb263044d3b435b479eb76ac674'])
        if len(target_list) < 1:
        #targets = set(['http://pathwaycommons.org/pc12/Protein_83baeb7dc5ecdcda2877b086aebb603f'])
            targets = set(['http://pathwaycommons.org/pc12/Complex_2f0a923dd4cf0b828c8176a962e14011'])
        else:
            targets = set(target_list)
    elif len(target_list) > 0:
        sources = set()
        targets = set(target_list)
    print(('sources,targets',sources,targets))
    ## backward star sources
    high_penalty_sources = set([n for n in H.get_node_set() if len(H.get_backward_star(n))==0]).difference(sources)
    #all_sources = set(sources.union(high_penalty_sources))
    all_sources = sources 
    if len(all_sources.intersection(targets)) >0:
        print('Warning: removing %d nodes from targets that were both sources and targets' % (len(set(sources).intersection(set(targets)))))
        targets = [t for t in targets if t not in all_sources]

    print('%d sources and %d targets; %d high penalty sources (empty backward stars)'  % (len(sources),len(targets),len(high_penalty_sources)))
    name_dict = {}
    return sources,targets,high_penalty_sources,name_dict


def add_super_nodes(H,sources,targets,high_penalty_sources,dataset,stoichiometries=None,single_targets_edge=False):

    super_source = 'SUPERSOURCE'
    print('-------------------------------sourcelen is: ',len(sources))
    largenum = len(H.get_hyperedge_id_set())+100
    #for s in sources:
        #H.add_hyperedge(set([super_source]), set([s]), weight=1)
    hid = H.add_hyperedge(set([super_source]),set(sources.union(high_penalty_sources)),weight=0)
    
    super_target = 'SUPERTARGET'
    for t in targets:
        hid = H.add_hyperedge(set([t]),set([super_target]), weight=0)
    return super_source,super_target

def median(vals):
    vals.sort()
    if len(vals) % 2 != 0:
        return vals[int(len(vals)/2)]
    else:
        return (vals[len(vals)/2] + vals[len(vals) / 2 + 1])/2

def parseOptions(args):
    desc = 'python master-script.py [options]'
    parser = OptionParser(usage=desc)

    # General Options
    parser.add_option('','--force',action='store_true',help='Overwrite files if they exist.')
    parser.add_option('','--printonly',action='store_true',help='Print commands to screen, but do not execute them.')
    
    # EXPERIMENTS/TESTS
    parser.add_option('','--name',type='string',default='WNT5A',help='Name of dataset (WNT5A, CTNNB1, WNT, or ALL). Default=WNT.')
    parser.add_option('','--cyclicparallel',action='store_true',help='Compute shortest cyclic hyperpath from S to T.')
    parser.add_option('','--ilpchar',type='string',default='7',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--pathwayreconstruct',action='store_true',help='enumerate heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--source',type='string',action='append',help='Sources. Default = WNT5A')
    parser.add_option('','--target',type='string',action='append',help='Targets. Default = CTNNB1')
    parser.add_option('','--pathway',type='string',action='append',help='Targets. Default = CTNNB1')


    opts,args = parser.parse_args()

    print('OPTIONS ARE',opts)
    return opts

 
if __name__ == '__main__':
    main(sys.argv)

