#!/usr/bin/python

## ISMB poster submission
import os
import os.path
import sys
from optparse import OptionParser
import time
import pickle as pkl
import traceback

# for plotting

from ilp import *
from pathheuristic import *
from cutfinder import *

from halp.directed_hypergraph import DirectedHypergraph
from halp.utilities import directed_graph_transformations
from halp.algorithms.directed_paths import b_visit

## http://manual.graphspace.org/en/latest/Programmers_Guide.html#graphspace-python
## http://manual.graphspace.org/projects/graphspace-python/en/latest/
#from graphspace_python.graphs.classes.gsgraph import GSGraph
#from graphspace_python.api.client import GraphSpace
#import graphspace_interface as interface

user = 'skrieger@email.arizona.edu'
password = 'ducktail'
#graphspace = GraphSpace('skrieger@email.arizona.edu', 'ducktail')
    

## GLOBAL VARIABLES
DELIM=';'
FORCE = False
PRINTONLY = False
ROOTDIR = '/groups/kece/skrieger/hypergraphs/biopax-hypergraph/biopax-hypergraph/test'

COLORS = {'blue':'#6F95EB',
        'red':'#EB856F',
        'orange':'#EBCD6F',
        'gray':'#C6C6C6',
        'black':'#000000',
        'white':'#FFFFFF'
        }

def main(args):
    opts = parseOptions(args)

    # get dataset
    #H,G_complex,G = getHypergraphs(opts.name)
    #if opts.type == 'hypergraph':
        #pass
    #elif opts.type == 'graph-with-complexes':
        #H = G_complex
    #elif opts.type == 'graph':
        #H = G
    #else:
        #print 'ERROR: graph type is not allowed.'
        #sys.exit()

    #H,_,__ = make_hypergraph(ROOTDIR+'/parsed/WNT5A')
    if opts.pathwayreconstruct:
        #read through hyperedges file picking out those corresponding to the ones from a pathway
        pathwayedges = []

        import requests
        pathwayrequest = requests.get('https://www.pathwaycommons.org/pc2/search.json?q=pathway:{}&datasource=reactome&type=conversion'.format(opts.pathway[0]))
        if 'json' in pathwayrequest.headers['Content-Type']:
            jsonobject = pathwayrequest.json()
            for edge in jsonobject['searchHit']:
                pathwayedges.append(edge['uri'])

        Hpathway,_,__ = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),select_edges=pathwayedges)
        pathwayedges = Hpathway.get_hyperedge_id_set()
        sources = set([v for v in Hpathway.node_iterator() if len(Hpathway.get_backward_star(v)) == 0])
        targets = set([v for v in Hpathway.node_iterator() if len(Hpathway.get_forward_star(v)) == 0])
        node_dict = pkl.load(open('node_dict.pkl','rb'))
        print(len(sources),[node_dict[s] for s in sources])
        print(len(targets),[node_dict[t] for t in targets])
        print('number of hyperedges in the reactome pathway: ',len(pathwayedges))
        #if opts.pathway[0] ==

    if not opts.stoichiometry and not opts.sbml:
        beforehypergraphtime = time.time()
        H,_,__ = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name))
        afterbuildtime = time.time()
        print('time to build hypergraph', afterbuildtime - beforehypergraphtime)
    #for e in H.hyperedge_id_iterator():
        #print(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e))
    ## add a super souce and a super sink
        if not opts.pathwayreconstruct:
            print(('target: {0}'.format(opts.target)))
            a = H.get_backward_star(opts.target[0])
            H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
            source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name)
        elif opts.pathwayreconstruct:
            H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',sources,targets,keep_sources=True)
            source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,single_targets_edge=True)
        afteraltertime = time.time()
        print('time to alter hypergraph', afteraltertime - afterbuildtime)
    #print 'this is the source and target',H_sources,H_targets
    #source = opts.source

    # run
    if opts.acyclic:
        #added by Spencer
        H, vertex_dict = convert_hypergraph_nodes(H)
        #done added by Spencer
        ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #added by Spencer
        #tfile = open('targets.pkl','rb')
        #targetlist = pkl.load(tfile)
        #objectives = []
        #for t in targetlist:
            #H.remove_node('SUPERTARGET')
            #H.add_node('SUPERTARGET',{'label': 'SUPERTARGET'})
            #H.add_hyperedge([vertex_dict[t.strip()]],['SUPERTARGET'])
            ##this needs to be converted 
        ##done added by Spencer
        if opts.force or not os.path.isfile(ilpname):
            print(target,'target')
            make_shortest_acyclic_hyperpath_ilp(H,source,target,ilpname)
        else:
            print('not writing ILP. Use --force to overwrite.')
        variables,objective = runILP_singlesol(H,ilpname,outprefix,opts.force,target,source)
        writeObjective(H,open('objectivevars.ilp','w'))
        pathedges = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
        for v in variables[0]:
            if variables[0][v] > 0.2 and 'e' in v:
                print(v,variables[0][v])
        #objectives.append((objective,targetlist[t][0]))
        #print(objective,targetlist[t][0])
        #objfile = open('objectivepairsncipid.pkl','wb')
        #pkl.dump(objectives,objfile)
        if opts.viz:
            #vizHypergraph(H,'%s-%s-acyclic' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict)
            vizHypergraph(H,'%s-%s-acyclic-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True)

    if opts.oldcyclic:
        ilpname = '%s/results/reactome-%s-%s-cyclic.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-cyclic' % (ROOTDIR,opts.name,opts.type)
        if opts.force or not os.path.isfile(ilpname):
            make_shortest_cyclic_hyperpath_ilp(H,source,target,ilpname)
        else:
            print('not writing ILP. Use --force to overwrite.')
        variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target,source)
        if opts.viz:
            #vizHypergraph(H,'%s-%s-acyclic' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict)
            vizHypergraph(H,'%s-%s-acyclic-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True)

    if opts.cyclicparallel or opts.pathwayreconstruct:
        ilpname = '%s/results/reactome-%s-%s-cyclic.lp%s' % (ROOTDIR,opts.name,opts.type,opts.ilpchar)
        node_dict = pkl.load(open('node_dict.pkl','rb'))
        '''
        nodes = ['http://pathwaycommons.org/pc12/Complex_7eab16f802b2b126d178667fee01efb6',
                'http://pathwaycommons.org/pc12/Complex_2dd50360b547faa4b5b28bcefee8e8f2',
                'http://pathwaycommons.org/pc12/Complex_c71c1c1f69bf4b35779a8f2ad3843773',
                'http://pathwaycommons.org/pc12/Complex_f591639e20599a930a088b836e347b67',
                'http://pathwaycommons.org/pc12/Protein_3a0883385f11b4fd86bf783b66cc81d8',
                'http://pathwaycommons.org/pc12/Protein_7e62be6a609831a84e78885abee209e1']
        nodes = ['http://pathwaycommons.org/pc12/SmallMolecule_8fa502e6d896a6f4b4b89e6945f5293c',
                'http://pathwaycommons.org/pc12/SmallMolecule_842865be194ec2fa337762fcd4fa0df8',
                'http://pathwaycommons.org/pc12/SmallMolecule_706e9054935f8ed25566e2642e14aa0b',
                'http://pathwaycommons.org/pc12/SmallMolecule_9099c7ae4555aad6bee725f49b31bc8d',
                'http://pathwaycommons.org/pc12/SmallMolecule_1c32f39a947e8dc821e215cc0ffb819c']
        nodes = ['http://pathwaycommons.org/pc12/SmallMolecule_10473b9e6fef48ad318fe1f6089fe5b9']
        nodes = ['http://pathwaycommons.org/pc12/Complex_a8d85da8ff183e6ea3ac1e74ff81aa87',
        'http://pathwaycommons.org/pc12/SmallMolecule_7f377e4f772d39a140e59200ddda2b0d',
        'http://pathwaycommons.org/pc12/SmallMolecule_1eef4b089daf104d1c8d5fbe4c04970e',
        'http://pathwaycommons.org/pc12/SmallMolecule_546b0586c4fbf47215c2f825b02b24c9',
        'http://pathwaycommons.org/pc12/SmallMolecule_6c7327deaebe9f6883d665bc7ea8d62e',
        'http://pathwaycommons.org/pc12/SmallMolecule_2204e7c2401d0ab559e09bc62a178e29',
        'http://pathwaycommons.org/pc12/SmallMolecule_8d6aa72ba00b4e24cc9b0dc42bf6e0f9',
        'http://pathwaycommons.org/pc12/SmallMolecule_30e0cfd9fe691d042ae9b34c56e94a72',
        'http://pathwaycommons.org/pc12/SmallMolecule_c5a48e3e444e9dc3b2c3923cce0bfe90',
        'http://pathwaycommons.org/pc12/SmallMolecule_1c77b26a91a24e64e9249a5d202281a3',
        'http://pathwaycommons.org/pc12/SmallMolecule_5c7e3625bcca8ef4e7e1e09932cd0e5f',
        'http://pathwaycommons.org/pc12/SmallMolecule_3cf16238f3a29742331da93a6397eab8',
        'http://pathwaycommons.org/pc12/SmallMolecule_4cb51309f77ffb5ae10a04868937e734',
        'http://pathwaycommons.org/pc12/SmallMolecule_97506ef167dacf19b5d35959e0d1ae4a',
        'http://pathwaycommons.org/pc12/SmallMolecule_d9b61b35e3b4790a37a93805c8a5d6ff',
        'http://pathwaycommons.org/pc12/SmallMolecule_b26d290c977d5862092cd111e78343d7',
        'http://pathwaycommons.org/pc12/SmallMolecule_aca6dc292d4ebd7a30a9b3ebced5b7c0',
        'http://pathwaycommons.org/pc12/SmallMolecule_30f83eba6c184bd9a90626def02375df',
        'http://pathwaycommons.org/pc12/SmallMolecule_ea368c90d5a68cdfa5846928f1f714c1',
        'http://pathwaycommons.org/pc12/Complex_8635539c551e9f5c2e90b35ff40e221b',
        'http://pathwaycommons.org/pc12/SmallMolecule_5b7ebe5b70cc98f4780ea6e27afb42d5']
        nodes = ['http://pathwaycommons.org/pc12/SmallMolecule_5c7e3625bcca8ef4e7e1e09932cd0e5f',
                'http://pathwaycommons.org/pc12/Complex_a8d85da8ff183e6ea3ac1e74ff81aa87',
                'http://pathwaycommons.org/pc12/SmallMolecule_3cf16238f3a29742331da93a6397eab8',
        'http://pathwaycommons.org/pc12/Complex_8635539c551e9f5c2e90b35ff40e221b']
        for n in nodes:
            print(node_dict[n])
        '''
        H2 = get_doubly_reachable_graph(H,source,target,node_dict)
        print(len(H2.get_hyperedge_id_set()))
        sourcecut = set([source])
        sinkcut = H2.get_node_set()
        if target not in H2.get_node_set():
            print('target not reachable from source')
        sinkcut.remove(target)
        vertexheadcuts = (sourcecut,sinkcut)
        edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H2,vertexheadcuts,node_dict)
        #TODO: remove this
        #for e in pathwayedges:
            #print(e,H2.get_hyperedge_id(Hpathway.get_hyperedge_tail(e),Hpathway.get_hyperedge_head(e)))

        if opts.pathwayreconstruct:
            taildistancelist,Hnew,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
            #Hnew,P = find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,nodb=opts.nodb,notchh=opts.notchh,notc=opts.notc,nohh=opts.nohh,notarget=opts.notargetconstraint,noaugment=opts.noaugment,returnP=True)
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
            find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,nodb=opts.nodb,notchh=opts.notchh,notc=opts.notc,nohh=opts.nohh,notarget=opts.notargetconstraint,noaugment=opts.noaugment)
        #if opts.nodb:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,nodb=True)
        #elif opts.notchh:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,notchh=True)
        #elif opts.notc:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,notc=True)
        #elif opts.nohh:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,nohh=True)
        #elif opts.notargetconstraint:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,notarget=True)
        #elif opts.noaugment:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts,noaugment=True)
        #else:
            #find_shortest_hyperpath_parallel(H2,H2.get_node_set(),ilpname,target,source,node_dict,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,vertexheadcuts=vertexheadcuts)


    if opts.cyclic:
        ilpname = '%s/results/reactome-%s-%s-cyclic.lp%s' % (ROOTDIR,opts.name,opts.type,opts.ilpchar)
        if opts.force or not os.path.isfile(ilpname):
            #node_dict = printremappededges(H)
            node_dict = pkl.load(open('node_dict.pkl','rb'))
            afternodedictloadtime = time.time()
            print('time to load node dict', afternodedictloadtime - afteraltertime)
            print('running heuristic')
            taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
            print('done running heuristic')
            print('heuristic pathlength is {}'.format(weight(H2,P)))
            #print('targetback: ', H2.get_backward_star(target))
            #print('Path: ',P)
            print('getting cuts')
            edgeheadcuts, vertexheadcuts = find_cuts(H2,taildistancelist,node_dict,tailpathlist,head_cuts_only=True)
            print('done getting cuts')
            #pkl.dump(edgeheadcuts,open('edgeheadcuts.pkl','wb'))
            #pkl.dump(vertexheadcuts,open('vertexheadcuts.pkl','wb'))
            #pkl.dump(H2,open('H2.pkl','wb'))
            #pkl.dump(P,open('P.pkl','wb'))

            #target = 'SUPERTARGET'
            #source = 'SUPERSOURCE'

            #edgeheadcuts = pkl.load(open('edgeheadcuts.pkl','rb'))
            #vertexheadcuts = pkl.load(open('vertexheadcuts.pkl','rb'))
            #H2 = pkl.load(open('H2.pkl','rb'))
            #P = pkl.load(open('P.pkl','rb'))
            if opts.newinequalities:
                find_shortest_hyperpath_with_new_inequalities(H2,H2.get_node_set(),ilpname,target,source,headcuts=edgeheadcuts,ilpchar=opts.ilpchar,heuristicpath=P,vertexheadcuts=vertexheadcuts)
            else:
                if opts.realvaluedcuts:
                    find_shortest_hyperpath_with_factory(H2,H2.get_node_set(),ilpname,target,source,headcuts=headcuts,ilpchar=opts.ilpchar,realvaluedcuts=True)
                else:
                    find_shortest_hyperpath_with_factory(H2,H2.get_node_set(),ilpname,target,source,headcuts=headcuts,ilpchar=opts.ilpchar)
            #find_shortest_hyperpath_with_factory(H2,H2.get_node_set(),ilpname,target,source)
            #find_shortest_hyperpath_with_factory(H,H.get_node_set(),ilpname,target,source)
        else:
            print('not writing ILP. Use --force to overwrite.')
        #variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target,source)

    if opts.simplepath:
        ilpname = '%s/results/reactome-%s-%s-simplewalk.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-simplewalk' % (ROOTDIR,opts.name,opts.type)
        if opts.force or not os.path.isfile(ilpname):
            make_shortest_hyperpath_ilp_simple(H,source,target,ilpname)
        else:
            print('not writing ILP. Use --force to override.')
        variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)
        if opts.viz:
            #vizHypergraph(H,'%s-%s-simplewalk' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict)
            vizHypergraph(H,'%s-%s-simplewalk-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True)

    if opts.minimal_precursor:
        H,_,__,stoichiometries,negative_regulators = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,keep_negative_regulators=True)
        #print('numhyperedges',len(H.get_hyperedge_id_set()))
        #print(len(stoichiometries))
        #print(stoichiometries,negative_regulators)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        #source,target,stoichiometries = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,stoichiometries=stoichiometries)
        ilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        ilpfile = open(ilpname,'w')
        write_binary = True

        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved first precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first precursor ILP\n')
        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved precursor conservation MILP\n')
        except:
            print('not able to solve precursor conservation ILP\n')


        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved first negative precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first negative precursor ILP\n')
        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,negative_regulators=negative_regulators)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved negative precursor conservation MILP\n')
        except:
            print('not able to solve negative precursor conservation ILP\n')

    if opts.figure:

        H,_,__,stoichiometries,negative_regulators = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,keep_negative_regulators=True)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        if opts.figure2:
            pilpname = '%s/results/reactome-%s-%s-precursors2.lp' % (ROOTDIR,opts.name,opts.type)
            poutprefix = '%s/results/reactome-%s-%s-precursors2' % (ROOTDIR,opts.name,opts.type)

            milpname = '%s/results/reactome-%s-%s-min3.lp' % (ROOTDIR,opts.name,opts.type)
            moutprefix = '%s/results/reactome-%s-%s-min3' % (ROOTDIR,opts.name,opts.type)
        elif opts.figure3:
            pilpname = '%s/results/reactome-%s-%s-precursors%d.lp' % (ROOTDIR,opts.name,opts.type,opts.fignum)
            poutprefix = '%s/results/reactome-%s-%s-precursors%d' % (ROOTDIR,opts.name,opts.type,opts.fignum)

            milpname = '%s/results/reactome-%s-%s-min%d.lp' % (ROOTDIR,opts.name,opts.type,opts.fignum)
            moutprefix = '%s/results/reactome-%s-%s-min%d' % (ROOTDIR,opts.name,opts.type,opts.fignum)
        else:
            pilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
            poutprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)

            milpname = '%s/results/reactome-%s-%s-min2.lp' % (ROOTDIR,opts.name,opts.type)
            moutprefix = '%s/results/reactome-%s-%s-min2' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries)
        write_binary = True
        begintime = time.time()

        #pilpfile = open(pilpname,'w')
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        #pilpfile.close()
        #try:
            #mpsones,mpssourceones = solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
            #print('solved first precursor MILP\n')
            ##print('neg\n')
        #except:
            #print('not able to solve first precursor ILP\n')

        milpfile = open(milpname,'w')

        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            if opts.figure2:
                _, sfones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=4)
            elif opts.figure3:
                _, sfones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=opts.fignum)
            else:
                _, sfones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=1)
            print('solved first MILP\n')
            #print('neg\n')
        except:
            traceback.print_exc()
            print('not able to solve first ILP\n')

        milpfile = open(milpname,'w')
        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            if opts.figure2:
                _, negones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=5)
            elif opts.figure3:
                _, negones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=opts.fignum+1)
            else:
                _, negones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary,solnum=2)
            print('solved neg first MILP\n')
            #print('neg\n')
        except:
            traceback.print_exc()
            print('failedsolve first MILP\n')
        
        #if len(sfones) == len(negones):
            ##print('lens are the same, exiting')
            #exit()

        if opts.figure3:
            sonegregones,sourcesused,i = second_order_negative_regulation(milpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True,randnum=opts.fignum)
            print('solutions lengths: ',len(sfones),len(negones),len(sonegregones))
        elif opts.figure2:
            llll = 1
            #sonegregones,sourcesused,i = second_order_negative_regulation(milpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True,randnum=3)
        else:
            sonegregones,sourcesused,i = second_order_negative_regulation(milpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True)
            print('solutions lengths: ',len(sfones),len(negones),len(sonegregones))


        #source,target = add_super_nodes(H,sources,targets,high_penalty_sources,opts.name)
        #H2, vertex_dict = convert_hypergraph_nodes(H)
        #ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        #outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #make_shortest_acyclic_hyperpath_ilp(H2,source,target,silpname)

#
        #variables,objective = runILP_singlesol(H2,silpname,soutprefix,opts.force,target,source)
        #print('True objective is ', objective-1)
        #pathones = []
        #if objective > -1:
            #pathones = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
            #for v in variables[0]:
                #if variables[0][v] > 0.2 and 'e' in v:
                    #print(v,variables[0][v])
                    ##print(H.get_hyperedge_tail(v[2:]),H.get_hyperedge_head(v[2:]))
        #sptime = time.time()

        #now to just print out each of the hyperedges with vertices that are easy to read?
        if opts.figure2:
            vertexdict = {}
            print(opts.target[0])
            vertexdict[opts.target[0]] = 'v0'
            vertexdictnum = 1
            print('sfones')
            for e in sfones:
                vtaillist = H.get_hyperedge_tail(unb(e))
                vheadlist = H.get_hyperedge_head(unb(e))
                vtranstaillist = []
                for v in vtaillist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtranstaillist.append(vertexdict[v])
                vtransheadlist = []
                for v in vheadlist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtransheadlist.append(vertexdict[v])
                print(e,vtranstaillist,vtransheadlist)
                if unb(e) in negative_regulators:
                    print(negative_regulators[unb(e)])
                
            #print('mpsones')
            #for e in mpsones:
            print('negones')
            for e in negones:
                vtaillist = H.get_hyperedge_tail(unb(e))
                vheadlist = H.get_hyperedge_head(unb(e))
                vtranstaillist = []
                for v in vtaillist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtranstaillist.append(vertexdict[v])
                vtransheadlist = []
                for v in vheadlist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtransheadlist.append(vertexdict[v])
                print(e,vtranstaillist,vtransheadlist)
                if unb(e) in negative_regulators:
                    print(negative_regulators[unb(e)])

            '''
            print('sonegones')
            for e in sonegregones:
                vtaillist = H.get_hyperedge_tail(unb(e))
                vheadlist = H.get_hyperedge_head(unb(e))
                vtranstaillist = []
                for v in vtaillist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtranstaillist.append(vertexdict[v])
                vtransheadlist = []
                for v in vheadlist:
                    if v not in vertexdict:
                        vertexdict[v] = 'v{}'.format(vertexdictnum)
                        vertexdictnum += 1
                    vtransheadlist.append(vertexdict[v])
                print(e,vtranstaillist,vtransheadlist)
                if unb(e) in negative_regulators:
                    print(negative_regulators[unb(e)])

            '''
            node_dict = pkl.load(open('node_dict.pkl','rb'))
            #print('pathones')
            #for e in pathones:
                #vtaillist = H.get_hyperedge_tail(unb(e))
                #vheadlist = H.get_hyperedge_head(unb(e))
                #vtranstaillist = []
                #for v in vtaillist:
                    #if v not in vertexdict:
                        #vertexdict[v] = 'v{}'.format(vertexdictnum)
                        #vertexdictnum += 1
                    #vtranstaillist.append(vertexdict[v])
                #vtransheadlist = []
                #for v in vheadlist:
                    #if v not in vertexdict:
                        #vertexdict[v] = 'v{}'.format(vertexdictnum)
                        #vertexdictnum += 1
                    #vtransheadlist.append(vertexdict[v])
                #print(e,vtranstaillist,vtransheadlist)
                #print(e,[node_dict[w] for w in vtaillist],[node_dict[x] for x in vheadlist])
            reversedict = ['' for p in range(80)]
            for v in vertexdict:
                if int(vertexdict[v][1:]) < 80:
                    reversedict[int(vertexdict[v][1:])] = (node_dict[v],v)
            for w in range(80):
                print(w,reversedict[w])
                    #print(vertexdict[v],v)
                    #print(vertexdict[v],node_dict[v])

    if opts.hybrid:
        H,_,__,stoichiometries,reversible_reactions = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,return_reversible=True)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        pilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        poutprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)

        milpname = '%s/results/reactome-%s-%s-min.lp' % (ROOTDIR,opts.name,opts.type)
        moutprefix = '%s/results/reactome-%s-%s-min' % (ROOTDIR,opts.name,opts.type)

        silpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        soutprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries)
        write_binary = True
        begintime = time.time()
        ones = []

        milpfile = open(milpname,'w')

        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,sourcevars=True)
        milpfile.close()
        try:
            __, ones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved first MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first ILP\n')
        sftime1 = time.time()
        print('time for first sf MILP {}'.format(sftime1-begintime))

        if ones != []:

            pilpfile = open(pilpname,'w')
            write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions,edgenum=len(ones))
            #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
            pilpfile.close()
            sourceones = []
            try:
                _, sourceones = solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
                print('solved first precursor MILP\n')
                #print('neg\n')
            except:
                print('not able to solve first precursor ILP\n')
            minpretime1 = time.time()
            print('time for first mps MILP {}'.format(minpretime1-sftime1))
        else:
            print('did not attempt to solve mpsMILP')

            

    if opts.sbml:
        H,_,__,stoichiometries,reversible_reactions = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,return_reversible=True)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        pilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        poutprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)

        milpname = '%s/results/reactome-%s-%s-min.lp' % (ROOTDIR,opts.name,opts.type)
        moutprefix = '%s/results/reactome-%s-%s-min' % (ROOTDIR,opts.name,opts.type)

        silpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        soutprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries)
        write_binary = True
        begintime = time.time()

        pilpfile = open(pilpname,'w')
        write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions)
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        pilpfile.close()
        try:
            solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
            print('solved first precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first precursor ILP\n')
        minpretime1 = time.time()
        print('time for first precursor MILP {}'.format(minpretime1-begintime))
        pilpfile = open(pilpname,'w')
        write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=reversible_reactions)
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=False)
        pilpfile.close()
        try:
            solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
            print('solved precursor conservation MILP\n')
        except:
            print('not able to solve precursor conservation ILP\n')
        minpretime2 = time.time()
        print('time for conservation precursor MILP {}'.format(minpretime2-minpretime1))
        milpfile = open(milpname,'w')

        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved first MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first ILP\n')
        sftime1 = time.time()
        print('time for first sf MILP {}'.format(sftime1-minpretime2))

        milpfile = open(milpname,'w')
        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=reversible_reactions)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved conservation MILP\n')
        except:
            print('not able to solve conservation ILP\n')
        sftime2 = time.time()
        print('time for conservation sf MILP {}'.format(sftime2-sftime1))



        source,target = add_super_nodes(H,sources,targets,high_penalty_sources,opts.name)
        H, vertex_dict = convert_hypergraph_nodes(H)
        ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        make_shortest_acyclic_hyperpath_ilp(H,source,target,silpname)


        variables,objective = runILP_singlesol(H,silpname,soutprefix,opts.force,target,source)
        print('True objective is ', objective-1)
        if objective > -1:
            pathedges = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
            for v in variables[0]:
                if variables[0][v] > 0.2 and 'e' in v:
                    print(v,variables[0][v])
                    #print(H.get_hyperedge_tail(v[2:]),H.get_hyperedge_head(v[2:]))
        sptime = time.time()
        print('time for shortest path MILP {}'.format(sptime-sftime2))

    if opts.stoichiometry:
        timebegin = time.time()
        H,_,__,stoichiometries,negative_regulators = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,keep_negative_regulators=True)
        #print('numhyperedges',len(H.get_hyperedge_id_set()))
        #print(len(stoichiometries))
        #print(stoichiometries,negative_regulators)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        #source,target,stoichiometries = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,stoichiometries=stoichiometries)
        ilpname = '%s/results/reactome-%s-%s-stoichiometry2.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-stoichiometry2' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        if opts.second_order_neg_reg:
            timebegin = time.time()
            solution,sourcesused,i = second_order_negative_regulation(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            if solution != '':
                print('solved second order neg reg')
            else:
                print('failed to find second order neg reg solution')

            timeaccend = time.time()
            print('time for second order negative regulation {}'.format(timeaccend-timebegin))

            solution,sourcesused,i = second_order_negative_regulation(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=False)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            if solution != '':
                print('solved second order neg reg conservation')
            else:
                print('failed to find second order neg reg solution conservation')

            timeconsend = time.time()
            print('time for second order negative regulation {}'.format(timeconsend-timeaccend))

        elif opts.facenumeration:

            timebegin = time.time()
            solution,sourcesused,i = enumerate_minimal_factories(ilpname + '88',H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            timeend = time.time()
            print('time taken is {}'.format(timeend-timebegin))
            print('__________________ starting negative regulation 1st order __________________')
            solution,sourcesused,i = enumerate_minimal_factories(ilpname + '88',H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True,negative=True)

        elif opts.checknegreg:
            #Here we want to see if all of the factories and stuff are obeying 1st order negative regulation
            hyperpaths = parsehyperpaths(open('resultsstoichiometryheur.txt','r').readlines())
            minsourcefacs = parseminsourcefacs(open('resultsstoichiometryminp.txt','r').readlines())
            minedgefacs = parseminedgefacs(open('resultsstoichiometryedgeneg.txt','r').readlines())

            #now that we have the pathways, we need to determine how many of them violate neg reg
            validpcount = 0
            for p in hyperpaths:
                if validatepathwaynegreg(H,p,negative_regulators) == True:
                    validpcount += 1
            
            validmscount = 0
            for p in minsourcefacs:
                if validatepathwaynegreg(H,p,negative_regulators) == True:
                    validmscount += 1

            validmecount = 0
            for p in minedgefacs:
                if validatepathwaynegreg(H,p,negative_regulators) == True:
                    validmecount += 1

            print('valid hyperpath count: {}\n percent: {},{}'.format(validpcount,validpcount/len(hyperpaths),len(hyperpaths)))
            print('valid minsource count: {}\n percent: {},{}'.format(validmscount,validmscount/len(minsourcefacs),len(minsourcefacs)))
            print('valid minedge count: {}\n percent: {},{}'.format(validmecount,validmecount/len(minedgefacs),len(minedgefacs)))
        elif opts.checksymdiff:
            minsourcefacs0th = parseminsourcefacs(open('resultsstoichiometryminp.txt','r').readlines(),returndict=True)
            minedgefacs0th = parseminedgefacs(open('resultsstoichiometryedgeneg.txt','r').readlines(),returndict=True)
            minsourcefacs1st = parseminsourcefacs(open('resultsstoichiometryminp.txt','r').readlines(),returndict=True,negative=True)
            minedgefacs1st = parseminedgefacs(open('resultsstoichiometryedgeneg.txt','r').readlines(),returndict=True,negative=True)
            #Now just need the symmetric difference finder, which shouldn't be hard at all!
            edgesymdiffs = []
            sourcesymdiffs = []
            epercents = []
            spercents = []
            for t in minsourcefacs0th:
                if t in minsourcefacs1st:
                    sfacset = set(minsourcefacs0th[t])
                    sfacset1 = set(minsourcefacs1st[t])
                    sourcesymdiffs.append(len(sfacset.symmetric_difference(sfacset1)))
                    spercents.append(len(sfacset.intersection(sfacset1))/ max(len(sfacset),len(sfacset1)))
                    if len(sfacset.symmetric_difference(sfacset1)) > 40:
                        print(sfacset,sfacset1,len(sfacset),len(sfacset1))
            for t in minedgefacs0th:
                if t in minedgefacs1st:
                    efacset = set(minedgefacs0th[t])
                    efacset1 = set(minedgefacs1st[t])
                    edgesymdiffs.append(len(efacset.symmetric_difference(efacset1)))
                    epercents.append(len(efacset.intersection(efacset1))/ max(len(efacset),len(efacset1)))

            print('min source sym diff: {} max source sym diff: {} num nonzero sym diff: {},total len: {},min intersection/maxsize: {}'.format(min(sourcesymdiffs),max(sourcesymdiffs),sum([1 if sym > 0 else 0 for sym in sourcesymdiffs]),len(sourcesymdiffs),min(spercents)))
            print('min edge sym diff: {} max edge sym diff: {} num nonzero sym diff: {},total len: {},min intersection/maxsize: {}'.format(min(edgesymdiffs),max(edgesymdiffs),sum([1 if sym > 0 else 0 for sym in edgesymdiffs]),len(edgesymdiffs),min(epercents)))
        elif opts.symdiffoneinstance:
            for i in range(1,7):
                minedgefacs0th,minedgefacs1st = parseminedgefacsoneinstance(open('{}negenumresults.txt'.format(i),'r').readlines())
                print(len(minedgefacs0th),len(minedgefacs1st))
                #Now just need the symmetric difference finder, which shouldn't be hard at all!
                edgesymdiffs = []
                sourcesymdiffs = []
                epercents = []
                spercents = []
                maxsymdiff = 0
                maxinst = None
                for f0 in minedgefacs0th:
                    for f1 in minedgefacs1st:
                        efacset = set(f0)
                        efacset1 = set(f1)
                        symdiff = len(efacset.symmetric_difference(efacset1))
                        if symdiff > maxsymdiff:
                            maxsymdiff = symdiff
                            maxinst = [efacset,efacset1]
                        edgesymdiffs.append(symdiff)
                        epercents.append(len(efacset.intersection(efacset1))/ max(len(efacset),len(efacset1)))

                print('min edge sym diff: {} max edge sym diff: {} num nonzero sym diff: {},total len: {},min intersection/maxsize: {},num 0 order facs: {}, num 1 order facs: {}'.format(min(edgesymdiffs),max(edgesymdiffs),sum([1 if sym > 0 else 0 for sym in edgesymdiffs]),len(edgesymdiffs),min(epercents),len(minedgefacs0th),len(minedgefacs1st)))
                print(maxinst)
                print('\n\n')
        else:
            timebegin = time.time()
            ilpfile = open(ilpname,'w')
            write_binary = True
            #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
            write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
            ilpfile.close()
            try:
                solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                print('solved first MILP\n')
                #print('neg\n')
            except:
                
                print('not able to solve first ILP\n')
            ilpfile = open(ilpname,'w')
            write_binary = True
            #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary)
            write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,negative_regulators=negative_regulators)
            ilpfile.close()
            try:
                solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                print('solved conservation MILP\n')
            except:
                print('not able to solve conservation ILP\n')

            #source,target = add_super_nodes(H,sources,targets,high_penalty_sources,opts.name)
            #H, vertex_dict = convert_hypergraph_nodes(H)
            #done added by Spencer
            #ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
            #outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
            #make_shortest_acyclic_hyperpath_ilp(H,source,target,ilpname)


            #where time bound code goes
            #variables,objective = runILP_singlesol(H,ilpname,outprefix,opts.force,target,source)
            #print('True objective is ', objective-1)
            #if objective > -1:
                #pathedges = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
                #for v in variables[0]:
                    #if variables[0][v] > 0.2 and 'e' in v:
                        #print(v,variables[0][v])
                        #print(H.get_hyperedge_tail(v[2:]),H.get_hyperedge_head(v[2:]))
            
            #node_dict = pkl.load(open('node_dict.pkl','rb'))
            #taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
            #print(len(P),P)
            #make_shortest_hyperpath_ilp_simple(H,source,target,ilpname)
            #variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)
            


    if opts.heuristic:
        ilpname = '%s/results/reactome-%s-%s-heuristic.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-heuristic' % (ROOTDIR,opts.name,opts.type)
        #node_dict = printremappededges(H)
        #pkl.dump(node_dict,open('node_dict.pkl','wb'))
        node_dict = pkl.load(open('node_dict.pkl','rb'))
        #targetlist = [a.strip() for a in open('reactomecycles.txt','r').readlines()]
        test_source_sink_pairs_for_cycles(H,node_dict,opts,firsttimetargets=True)
        #test_source_sink_pairs_for_cycles(H,node_dict,opts,firsttimetargets=False)


        taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
        #taildistancelist,H2,P,recoveredlist,recoveredlistall = tail_path_heuristic(H,source,target,node_dict=node_dict,sinkrecover=True)
        #for edge in recoveredlistall:
            #if 6145 in H2.get_hyperedge_head(edge):
                #print(edge)
        #print('size recoverdlist',len(recoveredlist),len(recoveredlistall),len(H2.get_hyperedge_id_set()))
        #for edge in pathedges:
            #try:
                #H2edge = H2.get_hyperedge_id(H.get_hyperedge_tail(edge[2:]),H.get_hyperedge_head(edge[2:]))
                #if H2edge not in recoveredlist:
                    #print('edge not in recoveredlist',edge[2:],H.get_hyperedge_head(edge[2:]),H.get_hyperedge_tail(edge[2:]))
                    #if H2edge not in recoveredlistall:
                        #print('edge was also not in recoveredlistall')
                #else:
                    #print(edge, ' is in the list')
            #except:
                ##means the edge wasn't in H2
                #print(edge, ' was not in H2',H.get_hyperedge_tail(edge[2:]))
        #pkl.dump(taildistancelist,open('taildistancelist.pkl','wb'))
        #pkl.dump(tailpathlist,open('tailpathlist.pkl','wb'))
        #taildistancelist = pkl.load(open('taildistancelist.pkl','rb'))
        #tailpathlist = pkl.load(open('tailpathlist.pkl','rb'))
        #pkl.dump(H2,open('H2.pkl','wb'))
        #pkl.dump(P,open('P.pkl','wb'))
        #H2 = pkl.load(open('H2.pkl','rb'))
        #P = pkl.load(open('P.pkl','rb'))
        #taildistancelist2 = {}
        #H2 = DirectedHypergraph()
        #H2.add_nodes(H.get_node_set())
        #enum = 1
        #for edge in taildistancelist:
            #H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge),weight=H.get_hyperedge_weight(edge))
            #taildistancelist2['e{}'.format(enum)] = taildistancelist[edge]
            #enum += 1
        #cuts = find_cuts(H2,taildistancelist2)
        #headcuts,augmentedcuts = find_cuts(H2,taildistancelist,node_dict,tailpathlist)
        #pkl.dump(headcuts,open('headcuts.pkl','wb'))
        #pkl.dump(augmentedcuts,open('augmentedcuts.pkl','wb'))
        #headcuts = pkl.load(open('headcuts.pkl','rb'))
        #augmentedcuts = pkl.load(open('augmentedcuts.pkl','rb'))
        #print('{} total cuts'.format(len(cuts)))
        #print(cuts)

        if opts.force or not os.path.isfile(ilpname):
            make_shortest_cyclic_hyperpath_ilp(H2,source,target,ilpname)
        else:
            print('not writing ILP. Use --force to override.')

        #variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)
        #run_LP_heuristics(H2,H2.get_node_set(),ilpname,outprefix,1,target,source,headcuts,augmentedcuts,targetname=opts.target[0])
        if opts.viz:
            #vizHypergraph(H,'%s-%s-simplewalk' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict)
            vizHypergraph(H,'%s-%s-heuristic-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True)

    if opts.cheat:

        
        orig_hids = [hid for hid in H.hyperedge_id_iterator()]
        cheat_hids = set()
        for hid in orig_hids:
            tail = H.get_hyperedge_tail(hid)
            if len(tail) > 1:
                head = H.get_hyperedge_head(hid)
                for t in tail:
                    cheat_hids.add(H.add_hyperedge(set([t]),head,{'weight':1}))
        
        ilpname = '%s/results/reactome-%s-%s-cheat.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-cheat' % (ROOTDIR,opts.name,opts.type)
        if opts.force or not os.path.isfile(ilpname):
            #make_shortest_hyperpath_ilp_simple_sets(H,H_sources,target,ilpname)
            make_shortest_hyperpath_ilp_simple(H,source,target,ilpname)
        else:
            print('not writing ILP. Use --force to override.')
        variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)#,cheats=cheat_hids)
        if opts.viz:
            #vizHypergraph(H,'%s-%s-cheats' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,cheat_set = cheat_hids)
            #vizHypergraph(H,'%s-%s-cheats-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True)
            vizHypergraph(H,'%s-%s-cheats-solution-only' % (opts.name,opts.type),H_sources,H_targets,variables,H_name_dict,solonly=True,cheat_set = cheat_hids)

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

def find_alt_paths(P,H,source,target,node_dict={}):
    pathcount = 0
    for e in P:
        #remove e from H
        tail = H.get_hyperedge_tail(e)
        head = H.get_hyperedge_head(e)
        print((tail,head))
        H.remove_hyperedge(e)
        _,__,P2 = tail_path_heuristic(H,source,target,node_dict=node_dict)
        print(P2)
        if len(P2) > 0:
            pathcount += 1
        #do something with the path? Each path returned should be unique, so just count when P is empty

        H.add_hyperedge(tail,head)
    print(('this path had {0} out of {1} paths'.format(pathcount,len(P))))
        

def convert_hypergraph_nodes(H):
    H2 = DirectedHypergraph()
    edgecount = len(H.get_hyperedge_id_set())
    edgenums = sorted(int(e[1:]) for e in H.get_hyperedge_id_set())
    vertexdict = {}
    nodenum = 1
    for node in H.node_iterator():
        if node != 'SUPERTARGET' and node != 'SUPERSOURCE':
            vertexdict[node] = nodenum
            nodenum += 1
        else:
            vertexdict[node] = node
    for e in edgenums:
        if H.has_hyperedge_id('e{0}'.format(e)):
            tail = [vertexdict[v] for v in H.get_hyperedge_tail('e{0}'.format(e))]
            for node in tail:
                if type(node) is not int:
                    print(node)
            head = [vertexdict[v] for v in H.get_hyperedge_head('e{0}'.format(e))]
            for node in head:
                if type(node) is not int:
                    print(node)
            H2.add_hyperedge(set(tail),set(head),weight=H.get_hyperedge_weight('e{0}'.format(e)))
    for node in H2.node_iterator():
        if type(node) is not int:
            print(node)
    for edge in H2.hyperedge_id_iterator():
        for node in H2.get_hyperedge_tail(edge):
            if type(node) is not int:
                print(node)
        for node in H2.get_hyperedge_head(edge):
            if type(node) is not int:
                print(node)
    return H2,vertexdict

def test_source_sink_pairs_for_cycles(H,node_dict,opts,targetlist=None,firsttimetargets=True):
    fileappendage = 'ecoli'
    if targetlist == None:
        #If it is not the first time for us getting the targets, we want to go through the checkpointing file 'checkedfor...' and find the spot we last left off on.
        if firsttimetargets == False:
            incheckedf = open('checkedforcycles{0}.txt'.format(fileappendage),'r')
            incheckedfile = incheckedf.readlines()
            if len(incheckedfile) > 0:
                lasttarget = incheckedfile[-1].strip()
                print(('lasttarget',lasttarget))
            else:
                lasttarget = ''
            incheckedf.close()
            checkedfile = open('checkedforcycles{0}.txt'.format(fileappendage),'w')
            #targetedge = H.get_backward_star('SUPERTARGET')
            #H.remove_node('SUPERTARGET')
        #Then we want to get the target list we want to process, NOTE: This is done each time from scratch
            targets = []
            for v in H.get_node_set():
                if len(H.get_forward_star(v)) == 0 and len(H.get_backward_star(v)) > 0 and v != 'SUPERTARGET':
                    targets.append(v)
            targets = list(set(targets))
            oldcyclesfile = open('cycles{0}.txt'.format(fileappendage),'r')
            oldlines = oldcyclesfile.readlines()
            oldcyclesfile.close()
        else:
            targets = []
            for v in H.get_node_set():
                if len(H.get_forward_star(v)) == 0 and len(H.get_backward_star(v)) > 0:
                    targets.append(v)
            smallf = open('{0}targets.txt'.format(fileappendage),'w')
            for t in targets:
                smallf.write('{0}\n'.format(t))
            print(targets)
            checkedfile = open('checkedforcycles{0}.txt'.format(fileappendage),'w')
    targetfile = open('cycles{0}.txt'.format(fileappendage),'w')
    if targetlist == None:
        #Here we write back all of the old cycle lines to the file
        lasttargetsspot = 0
        if firsttimetargets == False:
            for line in oldlines:
                targetfile.write(line)
            if lasttarget != '':
                lasttargetsspot = targets.index(lasttarget)
    else:
        targets = targetlist
        lasttargetsspot = 0
    for target in targets[lasttargetsspot:]:
        if targetlist == None:
            checkedfile.write('{0}\n'.format(target))
        H.add_node('SUPERTARGET',{'label': 'SUPERTARGET'})
        H.add_hyperedge([target],['SUPERTARGET'])
        print((target,H.get_backward_star(target)))
        print((targets.index(target),len(targets)))
        if opts.bestinedgelist:
            taildistancelist,H2,P = tail_path_heuristic(H,'SUPERSOURCE','SUPERTARGET',node_dict=node_dict,best_in_edges=True)
        else:
            taildistancelist,H2,P = tail_path_heuristic(H,'SUPERSOURCE','SUPERTARGET',node_dict=node_dict)
        #possibly need to sort the edges first
        newP = []
        taillist = set(['SUPERSOURCE'])
        #print('sorting path')
        while len(newP) < len(P):
            #print(P,newP)
            for e in P:
                if e not in newP:
                    reached = True
                    #check if edge e has been reached by the other edges
                    for v in H2.get_hyperedge_tail(e):
                        #print(v)
                        if v not in taillist:
                            reached = False
                            break
                            #print('edge {} is not reached by {}'.format(e,taillist))
                    if reached == True:
                        newP.append(e)
                        for v in H2.get_hyperedge_head(e):
                            taillist.add(v)
                    
        #now is when we want to do the cyclefinder
        taillist = set()
        for e in newP:
            for v in H2.get_hyperedge_tail(e):
                taillist.add(v)
            for v in H2.get_hyperedge_head(e):
                if v in taillist:
                #if v in taillist and v not in H2.get_hyperedge_tail(e):
                    targetfile.write('target was {0}\ncycle with {1}\n'.format(target,v))
                    targetfile.write('path was \n')
                    for f in newP:
                        #targetfile.write('tail {}, head {}\n'.format([H2.get_node_attribute(w,'label') for w in H2.get_hyperedge_tail(f)],[H2.get_node_attribute(w,'label') for w in H2.get_hyperedge_head(f)]))
                        #targetfile.write('tail {0}, head {1}\n'.format([node_dict[w] for w in H2.get_hyperedge_tail(f)],[node_dict[w] for w in H2.get_hyperedge_head(f)]))
                        targetfile.write('tail {0}, head {1}\n'.format([w for w in H2.get_hyperedge_tail(f)],[w for w in H2.get_hyperedge_head(f)]))
                    vertexdict = {}
                    vertexcount = 0
                    vertexset = set()
                    for f in newP:
                        for v in H2.get_hyperedge_tail(f):
                            if v not in vertexset:
                                vertexdict[v] = 'v{0}'.format(vertexcount)
                                vertexcount += 1
                                vertexset.add(v)
                        for v in H2.get_hyperedge_head(f):
                            if v not in vertexset:
                                vertexdict[v] = 'v{0}'.format(vertexcount)
                                vertexcount += 1
                                vertexset.add(v)
                    for f in newP:
                        targetfile.write('tail {0} head {1}\n'.format([vertexdict[w] for w in H2.get_hyperedge_tail(f)],[vertexdict[w] for w in H2.get_hyperedge_head(f)]))
                    print('\n\n\n\n\nPOTENTIAL CYCLE\n\n\n\n\n')
        for edge in H.get_backward_star('SUPERTARGET'):
            H.remove_hyperedge(edge)
        #ilpname = '%s/results/reactome-%s-%s-heuristic.lp' % (ROOTDIR,opts.name,opts.type)
        #outprefix = '%s/results/reactome-%s-%s-heuristic' % (ROOTDIR,opts.name,opts.type)
        #cuts,normalcuts,headcuts = find_cuts(H2,taildistancelist,node_dict)
        #if opts.force or not os.path.isfile(ilpname):
            #make_shortest_cyclic_hyperpath_ilp(H2,'SUPERSOURCE','SUPERTARGET',ilpname)
        #else:
            #print 'not writing ILP. Use --force to override.'

        #variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)
        #run_LP_heuristics(H2,H2.get_node_set(),ilpname,outprefix,1,'SUPERTARGET','SUPERSOURCE',cuts,normalcuts,headcuts,targetname=opts.target[0])
            
        


def runILP_singlesol(H,ilpname,outprefix,force,target,source,cheats=None):
    objective=-1
    if force or not os.path.isfile('%s-1.variables' % (outprefix)):
        if not cheats:
            numsols,numoptobjective,allvars,objective = solveILP(H,H.get_node_set(),ilpname,outprefix,1,target,source,printobjective=True)
        else:
            numsols,numoptobjective,allvars = solveILP_cheats(H,H.get_node_set(),ilpname,outprefix,1,cheats)
        #allvars = allvars[0]
    else:
        print('not solving ILP. Use --force to override.')
        allvars = {}
        with open('%s-1.variables' % (outprefix)) as fin:
            for line in fin:
                if line[0] == '#': 
                    continue
                row = line.strip().split()
                if 'o_' in row[0]:
                    allvars[row[0]] = float(row[1])
                else:
                    allvars[row[0]] = int(row[2])
    return allvars,objective
#############################################

def execute(command,file_to_check,display_name):
    if not FORCE and os.path.isfile(file_to_check):
        print('Skipping %s: file %s exists. Use --force to overwrite.' % (display_name,file_to_check))
        return
    print(command)
    if not PRINTONLY:
        os.system(command)
    return

def checkdir(path):
    if not os.path.isdir(path):
        print('Making directory %s' % (path))
        os.makedirs(path)
    return


def make_hypergraph(file_prefix,delim=';',sep='\t',keep_singleton_nodes=False,stoichiometry=False,keep_negative_regulators=False,return_reversible=False,target=None,select_edges=[]):
    hypernodes = {}
    with open(file_prefix+'-hypernodes.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split(sep)
            ## TODO fix -- this happens when a complex contains other complexes
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
    if stoichiometry == True:
        stoichiometries = {}
    if keep_negative_regulators == True:
        negative_regulators = {}
    if return_reversible == True:
        reversible_reactions = {}
        prevhead = set()
        prevtail = set()
        previd = None

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
                if keep_negative_regulators == True:
                    if row[3] != 'None':
                        negative_regulators[hid] = row[3].split(delim)
                if stoichiometry == True:
                    stoichiometries[hid] = {}
                    if row[5] != 'None':
                        for molstoi in row[5].split(delim):
                            molecule = 'http:{}'.format(molstoi.split(':')[1])
                            if 'parsed' in file_prefix.split('/')[-1]:
                                stoichio = molstoi.split(':')[1]
                            else:
                                stoichio = molstoi.split(':')[2]
                            stoichiometries[hid][molecule] = float(stoichio)
                    for molecule in tail:
                        if molecule not in stoichiometries[hid]:
                            stoichiometries[hid][molecule] = 1.0
                    for molecule in head:
                        if molecule not in stoichiometries[hid]:
                            stoichiometries[hid][molecule] = 1.0
                if return_reversible == True:
                    if head == prevtail and tail == prevhead:
                        reversible_reactions[hid] = previd
                        reversible_reactions[previd] = hid
                    prevhead = head
                    prevtail = tail
                    previd = hid

    if 'WNT' in file_prefix.split('/') :
        print(('num vertices {0}'.format(len(H.get_node_set()))))
        print(('num edges {0}'.format(len(H.get_hyperedge_id_set()))))
        print(('largest head {0}'.format(max(headsizes))))
        print(('median head {0}'.format(median(headsizes))))
        print(('largest tail {0}'.format(max(tailsizes))))
        print(('median tail {0}'.format(median(tailsizes))))
        print(('num selfloops {0}'.format(len(selfloops))))
        selfset = set()
        for l in selfloops:
            selfset.update(l)
        for v in selfset:
            intersection = H.get_backward_star(v).intersection(H.get_forward_star(v))
            if len(intersection) == len(H.get_backward_star(v)):
                noinselfloops += 1
        print(('num no in selfloops {0}'.format(noinselfloops)))
        for v in H.get_node_set():
            outdegree.append(len(H.get_forward_star(v)))
            indegree.append(len(H.get_backward_star(v)))
            if len(H.get_forward_star(v)) == 0:
                numtargets += 1
            if len(H.get_backward_star(v)) == 0:
                numsources += 1
        print(('max outdegree {0}'.format(max(outdegree))))
        print(('max indegree {0}'.format(max(indegree))))
        print(('median indegree {0}'.format(median(indegree))))
        print(('median outdegree {0}'.format(median(outdegree))))
        print(('num targets {0}'.format(numtargets)))
        print(('num sources {0}'.format(numsources)))
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

    #print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' %
        #(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

    if stoichiometry == False:
        if keep_negative_regulators == False:
            if return_reversible == False:
                return H, identifier2id, id2identifier
            else:
                return H, identifier2id, id2identifier, reversible_reactions

        else:
            if return_reversible == False:
                return H, identifier2id, id2identifier, negative_regulators
            else:
                return H, identifier2id, id2identifier, negative_regulators, reversible_reactions
    else:
        if keep_negative_regulators == False:
            if return_reversible == False:
                return H, identifier2id, id2identifier, stoichiometries
            else:
                return H, identifier2id, id2identifier, stoichiometries, reversible_reactions
        else:
            if return_reversible == False:
                return H, identifier2id, id2identifier, stoichiometries, negative_regulators
            else:
                return H, identifier2id, id2identifier, stoichiometries, negative_regulators, reversible_reactions
        

def getHypergraphs(name):
    '''
    Get hypergraph and graph-with-complexes.
    '''
    H = DirectedHypergraph()
    H.read('%s/parsed/reactome-%s-hypergraph_remapped.txt' % (ROOTDIR,name))
    
    G_complex = DirectedHypergraph()
    G_complex.read('%s/parsed/reactome-%s-graph-with-complexes_remapped.txt' % (ROOTDIR,name))

    G = DirectedHypergraph()
    G.read('%s/parsed/reactome-%s-graph_remapped.txt' % (ROOTDIR,name))

    '''
    hnode_names = {}
    hnode_types = {}
    ## annotate with common names and hypernodes
    if os.path.isfile('%s/Data/%s/combined-hypergraph/%s-hypernodes.txt' % (ROOTDIR,name,name)):
        hnodes = {}
        with open('%s/Data/%s/combined-hypergraph/%s-hypernodes.txt' % (ROOTDIR,name,name)) as fin:
            for line in fin:
                if line[0] == '#':
                    continue # skip header
                row = line.strip().split()
                hnodes[row[0]] = row[1].split(DELIM)
                hnode_types[row[0]] = row[2]
                hnode_names[row[0]] = row[3]
        
        for hnode in H.node_iterator():
            H.add_node(hnode,{'elements':hnodes[hnode],'name':hnode_names[hnode],'type':hnode_types[hnode]})

    ## Convert to Graph with Complexes
    G = directed_graph_transformations.to_graph_decomposition(H)
    ## TODO: Convert to Graph
    '''
    return H,G_complex,G

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
    elif name == 'WNT5A':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee'])
        #sources = set(['http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #example 2 to break the loop
        #sources = set(['http://pathwaycommons.org/pc12/Protein_bdf326eedb65f7fe6a0259c5cb8c4ed4','http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #Trying to find cyclic

        targets = set(['http://pathwaycommons.org/pc12/Protein_6eae9e1fb8e906b20a3ebdf4485a4a3d']) #example 1
        #targets = set(['http://pathwaycommons.org/pc12/Protein_bc47c96d7c652d22f94260b30d5c8043']) #example 2
        #targets = set(['http://pathwaycommons.org/pc12/Complex_46f99a13ca1b39c8d93926b9b394c395'])# trying to find cyclic
    elif name == 'WNT':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Protein_3c9a7ce5eec5c6ec97c1a3008c3c9c99','http://pathwaycommons.org/pc12/Protein_1625818ba862a465e6bfe45c1a57c0ec'])
    elif name == 'allpid':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Complex_81ba3b0707b6c6abd477dd59315147f4'])
        if ':' in target_list[0]:
            targets = set([target_list[0]])
    elif name == 'allreactome' or name == 'allreactomestoichiometry':
        sources = set(['http://pathwaycommons.org/pc12/Complex_eb845eb263044d3b435b479eb76ac674'])
        if len(target_list) < 1:
        #targets = set(['http://pathwaycommons.org/pc12/Protein_83baeb7dc5ecdcda2877b086aebb603f'])
            targets = set(['http://pathwaycommons.org/pc12/Complex_2f0a923dd4cf0b828c8176a962e14011'])
        else:
            targets = set(target_list)
    elif name == 'WNTSTOI2':
        sources = set(['s'])
        targets = set(['t'])
    elif name == 'WNTSTOI':
        sources = set()
        targets = set([v for v in H.node_iterator() if len(H.get_forward_star(v)) == 0])
    elif len(target_list) > 0:
        sources = set()
        targets = set(target_list)
    else:
        sources = set()
        targets = set([min(H.get_node_set())])

    print(('sources,targets',sources,targets))
    #for v in targets:
        #for e in H.get_backward_star(v):
            #print((e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e)))
    #for e in H.hyperedge_id_iterator():
        #for v in H.get_hyperedge_head(e):
            #if v in targets:
                #print(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e))

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


def getSourcesTargets(name,H,graph_type,source_list,target_list):
    sources = set()
    targets = set()
    name_dict = {} # {id: name}
    print('name of file','%s/parsed/reactome-%s-%s_nodes_remapped.txt' % (ROOTDIR,name,graph_type))
    with open('%s/parsed/reactome-%s-%s_nodes_remapped.txt' % (ROOTDIR,name,graph_type)) as fin:
        for line in fin:
            row = line.strip().split('\t')
            name_dict[row[0]] = row[1]
            for s in source_list:
                if s.lower() in row[1].lower().split(':') and row[0] not in sources:
                #if s == row[1] and row[0] not in sources:
                    print('Source: ',s,row[1])
                    sources.add(row[0])
            for t in target_list:
                if t.lower() in row[1].lower().split(':')  and row[0] not in targets:
                #if t == row[1] and row[0] not in targets:
                    print('Target: ',t,row[1])
                    if row[1] != 'CTNNB1':
                        targets.add(row[0])
    print(sources,targets)

    ## add names (this will overwrite the attributes in the original names.)
    #for nid in H.node_iterator():
    #    if nid in name_dict:
    #        H.add_node(nid,{'name':name_dict[nid],'weight':1})

    ## backward star sources
    high_penalty_sources = set([n for n in H.get_node_set() if len(H.get_backward_star(n))==0]).difference(sources)
    #all_sources = set(sources.union(high_penalty_sources))
    all_sources = sources 
    if len(all_sources.intersection(targets)) >0:
        print('Warning: removing %d nodes from targets that were both sources and targets' % (len(set(sources).intersection(set(targets)))))
        targets = [t for t in targets if t not in all_sources]

    print('%d sources and %d targets; %d high penalty sources (empty backward stars)'  % (len(sources),len(targets),len(high_penalty_sources)))
    return sources,targets,high_penalty_sources,name_dict

def add_super_nodes(H,sources,targets,high_penalty_sources,dataset,stoichiometries=None,single_targets_edge=False):

    super_source = 'SUPERSOURCE'
    print('-------------------------------sourcelen is: ',len(sources))
    largenum = len(H.get_hyperedge_id_set())+100
    #for s in sources:
        #H.add_hyperedge(set([super_source]), set([s]), weight=1)
    hid = H.add_hyperedge(set([super_source]),set(sources.union(high_penalty_sources)),weight=0)
    if stoichiometries != None:
        stoichiometries[hid] = {}
        stoichiometries[hid][super_source] = 1.0
        for s in set(sources.union(high_penalty_sources)):
            stoichiometries[hid][s] = 1.0
    
    #for s in high_penalty_sources:
    ##this is where I changed the weight of high-penalty sources -Anna
        #H.add_hyperedge(set([super_source]), set([s]), weight=1)    
        #H.add_hyperedge(set([super_source]), set([s]), weight=100)    

    super_target = 'SUPERTARGET'
    if dataset == 'test' or dataset == 'reactome' or dataset == 'WNT' or single_targets_edge == True:
        hid = H.add_hyperedge(set(targets),set([super_target]), weight=1)
        if stoichiometries != None:
            stoichiometries[hid] = {}
            stoichiometries[hid][super_target] = 1.0
            for s in set(targets):
                stoichiometries[hid][s] = 1.0
    else: # dataset == 'ncipid'
        for t in targets:
            hid = H.add_hyperedge(set([t]),set([super_target]), weight=1)
            if stoichiometries != None:
                stoichiometries[hid] = {}
                stoichiometries[hid][super_target] = 1.0
                stoichiometries[hid][t] = 1.0

    ## TODO: Convert to Graph
    if stoichiometries == None:
        return super_source,super_target
    else:
        return super_source,super_target,stoichiometries

def vizHypergraph(H,graphid,sources,targets,variables,name_dict,solonly=False,cheat_set=[]):
    #here we need to put something that only visualizes the nodes based on the edges that are present, instead of the nodes.
    print(variables)
    ansnodes = {}
    for v in list(variables.keys()):
        ansnodes[v] = variables[v]
        for s in H.get_hyperedge_head(una(v)):
            if variables[v] == 1.0:
                ansnodes[a(s)] = 1.0
            else:
                if a(s) not in list(ansnodes.keys()):
                    ansnodes[a(s)] = 0.0
                
    ansnodes['a_SUPERSOURCE'] = 1.0
    variables = ansnodes
    
    print(variables)
    if solonly:
        print([v for v in list(variables.keys()) if 'a_' in v and variables[v] == 1])
        varnames = set(['_'.join(v.split('_')[1:]) for v in list(variables.keys()) if 'a_' in v and variables[v] == 1])
        numnodes = len(varnames.intersection(H.get_node_set()))
        numedges = len(varnames.intersection(H.get_hyperedge_id_set()))
        if numnodes > 10:
            print('Not posting hypergraph SOLUTION with %d nodes and %d hyperedges' % (numnodes,numedges))
        print('Posting hypergraph SOLUTION with %d nodes and %d hyperedges' % (numnodes,numedges))
    else:
        if directed_statistics.number_of_nodes(H) > 10:
            print('Not posting hypergraph with %d nodes and %d hyperedges' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H)))
            return
        print('Posting hypergraph with %d nodes and %d hyperedges' % (directed_statistics.number_of_nodes(H),directed_statistics.number_of_hyperedges(H)))

    #G = GSGraph()
    ## first add nodes
    for name in H.node_iterator():
        #if 'name' in H.get_node_attributes(name).keys():
            #if variables['a_'+name] == 1:
        #        label = H.get_node_attribute(name,'name')+'\n(%s %s)' % (H.get_node_attribute(name,'type'),name)
        #    else:
        #    label = H.get_node_attribute(name,'name')
        #else:
        label = name_dict.get(name,name)
        popup = 'Node "%s" (id: %s) <br>' % (label,name)

        if sources and name in sources:
            shape='triangle'
        elif targets and name in targets:
            shape='rectangle'
        else:
            shape = 'ellipse'
        
        if variables['a_'+name] == 1:
            color = COLORS['blue']
            width=5+10*len(label)
            height=60
        else:
            if solonly:
                continue
            color = COLORS['gray']
            width=20
            height=20
        print('Adding node',name,label)
        G.add_node(name,label=label,popup=popup)
        G.add_node_style(name,color=color,shape=shape,width=width,height=height)
        
    ## then add hyperedges
    for name in H.hyperedge_id_iterator():
        if solonly and variables['a_'+name] == 0:
            continue

        print('Hyperedge %s connecting %s and %s' % (name,';'.join(H.get_hyperedge_tail(name)),';'.join(H.get_hyperedge_head(name))))
        popup='Hyperedge "%s"<br>Weight %f<br>Tails: %s<br>Heads: %s' % (name,H.get_hyperedge_attribute(name,'weight'),\
            ';'.join(H.get_hyperedge_tail(name)),';'.join(H.get_hyperedge_head(name)))

        if variables['a_'+name] == 1:
            if H.get_hyperedge_attribute(name,'weight') > 100 : ## largenum
                color = COLORS['red']
            elif name in cheat_set:
                color = COLORS['orange']
            else:
                color = COLORS['blue']
            width = 3
        else:
            color = COLORS['black']
            width = 1

        if len(H.get_hyperedge_tail(name)) == 1 and len(H.get_hyperedge_head(name))==1:
            t = [i for i in H.get_hyperedge_tail(name)][0]
            h = [i for i in H.get_hyperedge_head(name)][0]
            G.add_edge(t,h,popup=popup)
            G.add_edge_style(t,h,directed=True,color=color,width=width)
        else:
            G.add_node(name,label='',popup=popup)
            G.add_node_style(name,color=COLORS['white'],border_color=COLORS['black'],shape='rectangle',height=10,width=10)
        
            for tail in H.get_hyperedge_tail(name):
                G.add_edge(tail,name,popup=popup)
                G.add_edge_style(tail,name,directed=False,color=color,width=width)
            for head in H.get_hyperedge_head(name):
                G.add_edge(name,head,popup=popup)
                G.add_edge_style(name,head,directed=True,color=color,width=width)

    print('Posting hypergraph with id "%s"' % (graphid))
    G2 = graphspace.get_graph(graphid)
    if G2: # graph exists; update it
        ret = graphspace.update_graph(graphid,graph=G)
    else: # post new graph (G2 is None)
        ret = graphspace.post_graph(G)
    #print ret
    print('done posting graph',graphid)


    #group = 'Hypershrubs'
    #print 'Sharing graph with group "%s"' % (group)
    #graphspace.shareGraph(graphid,user=user,password=password,group=group)

    return 
def median(vals):
    vals.sort()
    if len(vals) % 2 != 0:
        return vals[int(len(vals)/2)]
    else:
        return (vals[len(vals)/2] + vals[len(vals) / 2 + 1])/2

def parsehyperpaths(lines):
    #Just need to get all the acyclic hyperpaths out of the file, don't even need to know the target
    hyperpaths = []
    currp = []
    inhyperpath = False
    for i,line in enumerate(lines):
        if 'True objective is' in line:
            currp = []
            inhyperpath = True
        elif 'OPTIONS ARE' in line:
            inhyperpath = False
            if len(currp) > 0:
                hyperpaths.append(currp)
        elif inhyperpath == True:
            currp.append(line.split('\'')[1][2:])
    return hyperpaths

def parseminsourcefacs(lines,returndict=False,negative=False):
    if returndict == True:
        minsourcefacs = {}
    else:
        minsourcefacs = []

    target = ''
    phrase = 'solved first precursor MILP'
    if negative == True:
        phrase = 'solved first negative precursor MILP'
    for i,line in enumerate(lines):
        if 'target:' in line:
            target = line.split('\'')[1]
        elif phrase in line:
            minfacunrefined = lines[i-4].split(',')[:-3]
            minfac = [l.split('\'')[1][2:] for l in minfacunrefined]
            if returndict == True:
                minsourcefacs[target] = minfac
            else:
                minsourcefacs.append(minfac)

    return minsourcefacs

def parseminedgefacsoneinstance(lines):
    minedgefacs0 = []
    minedgefacs1 = []

    phrase = 'Total (root+'
    negatives = False
    for i,line in enumerate(lines):
        if phrase in line:
            minfacunrefined = lines[i+1].split(',')
            minfac = [l.split('\'')[1][2:] for l in minfacunrefined]
            if negatives == True:
                minedgefacs1.append(minfac)
            else:
                minedgefacs0.append(minfac)
        elif 'starting negative regulators' in line:
            negatives = True
        #elif 'iterations needed' in line:
            #if negatives =
        elif 'OPTIONS ARE' in line:
            #This means it was preempted and we should reset everything
            minedgefacs0 = []
            minedgefacs1 = []
            negatives = False


    return minedgefacs0,minedgefacs1

def parseminedgefacs(lines,returndict=False,negative=False):
    if returndict == True:
        minedgefacs = {}
    else:
        minedgefacs = []

    target = ''
    phrase = 'solved first MILP'
    if negative == True:
        phrase = 'solved neg first MILP'
    for i,line in enumerate(lines):
        if 'target:' in line:
            target = line.split('\'')[1]
        elif phrase in line:
            minfacunrefined = lines[i-4].split(',')[:-3]
            if ' 0' not in lines[i-3]:
                minfac = [l.split('\'')[1][2:] for l in minfacunrefined]
                if returndict == True:
                    minedgefacs[target] = minfac
                else:
                    minedgefacs.append(minfac)

    return minedgefacs

def validatepathwaynegreg(H,p,negregs):
    #negregs[e] gives the things that negreg e. So we need to see if any of the edges that makes that thing are in the pathway
    for e in p:
        if e in negregs:
            verticesnegrege = negregs[e]
            for v in verticesnegrege:
                if v in H.get_node_set():
                    for f in H.get_backward_star(v):
                        if f in p:
                            print(p,f,H.get_hyperedge_head(f),v)
                            return False
    return True

def parseOptions(args):
    desc = 'python master-script.py [options]'
    parser = OptionParser(usage=desc)

    # General Options
    parser.add_option('','--force',action='store_true',help='Overwrite files if they exist.')
    parser.add_option('','--printonly',action='store_true',help='Print commands to screen, but do not execute them.')
    
    # EXPERIMENTS/TESTS
    parser.add_option('','--name',type='string',default='WNT5A',help='Name of dataset (WNT5A, CTNNB1, WNT, or ALL). Default=WNT.')
    parser.add_option('','--type',type='string',default='hypergraph',help='graph type: hypergraph, graph-with-complexes, or graph. Default=hypergraph.')
    parser.add_option('','--acyclic',action='store_true',help='Compute shortest acyclic hyperpath from S to T.')
    parser.add_option('','--cyclic',action='store_true',help='Compute shortest cyclic hyperpath from S to T.')
    parser.add_option('','--cyclicparallel',action='store_true',help='Compute shortest cyclic hyperpath from S to T.')
    parser.add_option('','--oldcyclic',action='store_true',help='Compute shortest cyclic hyperpath from S to T.')
    parser.add_option('','--stoichiometry',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--minimal_precursor',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--sbml',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--figure',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--figure2',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--figure3',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--checknegreg',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--checksymdiff',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--tchh',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--notchh',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--notc',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--nodb',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--nohh',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--notargetconstraint',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--noaugment',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--symdiffoneinstance',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--fignum',type='int',default=7,help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--ilpchar',type='string',default='7',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--second_order_neg_reg',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--heuristic',action='store_true',help='Compute heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--hybrid',action='store_true',help='Compute heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--newinequalities',action='store_true',help='Compute heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--enumeration',action='store_true',help='enumerate heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--facenumeration',action='store_true',help='enumerate heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--pathwayreconstruct',action='store_true',help='enumerate heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--simplepath',action='store_true',help='Compute shortest hyperpath from S to T.')
    parser.add_option('','--bestinedgelist',action='store_true',help='compute with just best inedge list')
    parser.add_option('','--realvaluedcuts',action='store_true',help='compute with just best inedge list')
    parser.add_option('','--cheat',action='store_true',help='add cheat hyperedges')
    parser.add_option('','--source',type='string',action='append',help='Sources. Default = WNT5A')
    parser.add_option('','--target',type='string',action='append',help='Targets. Default = CTNNB1')
    parser.add_option('','--pathway',type='string',action='append',help='Targets. Default = CTNNB1')

    # VISUALIZE
    parser.add_option('','--viz',action='store_true',help='Post to GraphSpace')

    opts,args = parser.parse_args()
       ## set FORCE and PRINTONLY variables
    global FORCE, PRINTONLY
    if opts.force:
       FORCE = True
    if opts.printonly:
       PRINTONLY = True

    if not opts.source:
       opts.source = ['WNT5A']
    if not opts.target:
       opts.target = ['CTNNB1']
    print('OPTIONS ARE',opts)
    return opts

 
if __name__ == '__main__':
    main(sys.argv)

