from xml.dom import minidom
import sys
import cplex
import itertools
import heapq
import time
import pickle as pkl

from ilp import *
from pathheuristic import weight
from pathheuristic import *

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit,visit
from halp.utilities import directed_graph_transformations

EPSILON=0.0001
CONSTANT=1000000

#######################################################################################################################################################################
'''
NEW CUTFINDING STUFF
Kept the other things in this file because they might be useful utility functions. 

'''

def find_cuts(H,taildistancelist,node_dict,tailpathlist,head_cuts_only=False):

    '''
    Finds the cuts based off the heuristic path's tail distance list.
    '''
    begin = time.time()
    #print(taildistancelist)
    #print('finding cuts')
    #normalcuts = find_normal_cuts(H,taildistancelist,node_dict) #normal cuts have the list of vertices in C and a distance value associated with them
    #normalcuttime = time.time()
    #print(('finding normal cuts took {}'.format(normalcuttime - begin)))
    #print(('{} normal cuts'.format(len(normalcuts))))

    headcuts = find_head_cuts(H,taildistancelist,node_dict)
    headcutstuples = set()
    for cut in headcuts:
        headcutstuples.add(tuple(cut))
    edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H,headcutstuples,node_dict)
    if head_cuts_only == True:
        return edgeheadcuts,headcutstuples
    headcuttime = time.time()
    print(('finding head cuts took {}'.format(headcuttime - begin)))
    print(('{} head cuts'.format(len(headcuts))))

    augmentedcuts = augment_head_cuts(H,headcuts,tailpathlist)
    pkl.dump(augmentedcuts,open('augmentedcuts2.pkl','wb'))
    #augmentedcuts = pkl.load(open('augmentedcuts2.pkl','rb'))
    print('augmentedlen',len(augmentedcuts))
    #print(augmentedcuts)
    flataugmentedcuts = []
    for c in augmentedcuts:
        if isinstance(c[0],list) or isinstance(c[0],tuple):
            for b in c:
                flataugmentedcuts.append(b)
        else:
            flataugmentedcuts.append(c)
            
    print('flataugmentedlen',len(flataugmentedcuts))
    #print(flataugmentedcuts)

    #finalcuts = set([tuple(n[0]) for n in normalcuts])
    #finalcuts.add(('SUPERSOURCE',))
    #edgenormalcuts = convert_vertex_cuts_to_edge_cuts(H,finalcuts,node_dict)
    edgeaugmentedcuts = convert_vertex_cuts_to_edge_cuts(H,flataugmentedcuts,node_dict)
    print('edgeflataugmentedlen',len(edgeaugmentedcuts))
    print('edgeheadlen',len(edgeheadcuts))
    edgeaugmentedcuts = [c for c in edgeaugmentedcuts if len(c) > 2]
    print([c for c in edgeaugmentedcuts if len(c) < 3])
    #edgefinalcuts = edgenormalcuts.union(edgeheadcuts)
    #return edgefinalcuts
    #print('normal cuts',edgefinalcuts)
    #pseudofinaledges = find_pseudofinal_edges(H)
    #print('pseudo-final edges',pseudofinaledges)
    #timetolerance = .01
    #print('found pseudo-final edges')
    #for pfe in pseudofinaledges:
        #now = time.time()
        #if now - begin > timetolerance:
            #break
        #smallestinedge = find_smallest_inedge(H,pfe,taildistancelist)
        #cutgapval = (taildistancelist[smallestinedge] + weight(H,[smallestinedge]),taildistancelist[pfe])
        #cutgapstart,cutgapend = find_cutgap_index(cutgapval,normalcuts)
        #for ci in range(cutgapstart,cutgapend):
            #nownow = time.time()
            #if nownow - begin > timetolerance:
                #break
            #for cj in range(0,ci):
                #insidenow = time.time()
                #if insidenow - begin > timetolerance:
                    #break
                ##need to make the cut augments unique
                #augmentedcuts = augment_cut(H,normalcuts,ci,cj,pfe)
                #headaugmentedcuts = augment_cut(H,headcuts,ci,cj,pfe)
                #edgeaugmentedcuts = convert_vertex_cuts_to_edge_cuts(H,augmentedcuts,node_dict)
                #headedgeaugmentedcuts = convert_vertex_cuts_to_edge_cuts(H,headaugmentedcuts,node_dict)
                ##print('pfe',pfe,'ci and cj',ci,cj,'pfe tail',H.get_hyperedge_tail(pfe))
                ##print('augmented',edgeaugmentedcuts)
                ##finalcuts = finalcuts.union(augmentedcuts)
                ##finalcuts = finalcuts.union(headaugmentedcuts)
                #edgefinalcuts = edgefinalcuts.union(edgeaugmentedcuts)
                #edgefinalcuts = edgefinalcuts.union(headedgeaugmentedcuts)
    ##print(len(finalcuts))
    ##edgecuts = convert_vertex_cuts_to_edge_cuts(H,finalcuts,node_dict)
    ##print('edgecuts',edgecuts)
    ##print('len difference',len(edgecuts),len(edgefinalcuts))
    #sourcecut = convert_vertex_cuts_to_edge_cuts(H,[['SUPERSOURCE']],node_dict)
    ##edgecuts = edgecuts.union(sourcecut)
    #edgefinalcuts = edgefinalcuts.union(sourcecut)

    end = time.time()
    print(('whole cut finding took {}'.format(end - begin)))

    return edgeheadcuts,edgeaugmentedcuts

def find_cuts_run_heuristic(H,source,target,node_dict,return_dict):

    '''
    Runs the heuristic and then Finds the cuts based off the heuristic path's tail distance list.
    '''

    taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True,doublyreachablesubgraph=True,print_paths=True)
    #pkl.dump(taildistancelist,open('tdlist.pkl','wb'))
    #pkl.dump(H2,open('H2pickle.pkl','wb'))
    #pkl.dump(P,open('Ppickle.pkl','wb'))
    #pkl.dump(tailpathlist,open('tplist.pkl','wb'))

    #taildistancelist = pkl.load(open('tdlist.pkl','rb'))
    #H2 = pkl.load(open('H2pickle.pkl','rb'))
    #P = pkl.load(open('Ppickle.pkl','rb'))
    #tailpathlist = pkl.load(open('tplist.pkl','rb'))


    headcuts = find_head_cuts(H2,taildistancelist,node_dict)
    headcutstuples = set()
    for cut in headcuts:
        headcutstuples.add(tuple(cut))
    #edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H2,headcutstuples,node_dict)
    i = 0
    for e in P:
        for l in e:
            return_dict[2][i] = l.encode()
            i += 1
    return_dict[0].value = weight(H2,P)
    #We now want to serialize the vertex cuts, since they cannot be reconstructed from the edge cuts
    #serializededgeheadcuts = serialize_headcuts(edgeheadcuts)
    serializedvertexheadcuts = serialize_headcuts(headcutstuples,vertexcuts=True)
    for i in range(len(serializedvertexheadcuts)):
        return_dict[1][i] = serializedvertexheadcuts[i].encode()
    #print(return_dict,'return_dict in cutfinder')
    return

def unserialize_path(serialP):

    P = []
    edge = []
    for l in serialP:
        if l == 'e':
            if edge != []:
                P.append(''.join(edge))
            edge = ['e']
        else:
            edge.append(l)
    return P

def serialize_headcuts(cuts,vertexcuts=False):
    #takes the list of tuples and turns it into a list of strings
    scuts = []
    for cut in cuts:
        scuts.append('(')
        for e in cut:
            for l in e:
                scuts.append(l)
            if vertexcuts == True:
                scuts.append(',')
        scuts.append(')')
    return scuts

def unserialize_vertex_headcuts(cuts):
    #takes the list of strings and turns it into a list of tuples, stopping at the uninhabited values
    hcuts = []
    cut = []
    vertex = []
    #print('serialcuts',cuts)
    for l in cuts:
        if l == '(':
            cut = []
        elif l == ',':
            if len(vertex) > 0:
                cut.append(''.join(vertex))
            vertex = []
        elif l == ')':
            if len(vertex) > 0:
                cut.append(''.join(vertex))
            hcuts.append(cut)
            #print('I appended a cut')
            cut = []
            vertex = []
        else:
            vertex.append(l)
    return hcuts

def unserialize_headcuts(cuts):
    #takes the list of strings and turns it into a list of tuples, stopping at the uninhabited values
    hcuts = []
    cut = []
    edge = []
    print('serialcuts',cuts)
    for l in cuts:
        if l == '(':
            cut = []
        elif l == 'e':
            if len(edge) > 0:
                cut.append(''.join(edge))
            edge = ['e']
        elif l == ')':
            if len(edge) > 0:
                cut.append(''.join(edge))
            hcuts.append(cut)
            cut = []
            edge = []
        else:
            edge.append(l)
    return hcuts

def edge_to_vertex_cuts(H,cuts):
    #Takes edge cuts and converts them into vertex cuts
    vertexcuts = []
    for cut in cuts:
        sourceside = set()
        for e in cut:
            for v in H.get_hyperedge_tail(e):
                sourceside.add(v)
        vertexcuts.append(list(sourceside))
    return vertexcuts


def convert_vertex_cuts_to_edge_cuts(H,cutlist):
    '''
    Takes cuts in a cutlist that are defined by vertices in H and converts them into cuts that are based on edges in H
    '''
    edgecuts = set()
    for cut in range(len(cutlist)):
        E = find_crossing_edges(H,cut,cutlist)
        edgecuts.add(tuple(E))
    return edgecuts
        

def augment_cut(H,normalcuts,ci,cj,pfe):
    '''
    For each vertex subset up to a certain size, this will split the edges that cross Cj into two subsets, F and F-complement.
    Then it will find the reachable set of vertices R that are reachable from F but not from F-complement
    It will then check if R has at least one tail vertex of pfe in it to continue
    Where it will put R in C-bar and add one cut that will be returned
    '''

    augmentedcuts = set()
    crossingedges = find_crossing_edges(H,cj,normalcuts)
    #print(crossingedges)
    edgesubsets = enumerate_edge_subsets(crossingedges)
    #print(len(edgesubsets))
    for edgesubset,complement in edgesubsets:
        R = find_reachable_vertices(H,edgesubset,complement)
        pfesourcevertices,doesintersect = check_pfe_intersection(H,pfe,R)
        if doesintersect == True:
            newcut = make_augmented_cut(normalcuts[cj][0],R,pfesourcevertices)
            augmentedcuts.add(newcut)

    return augmentedcuts


def make_augmented_cut(cut,R,C):
    '''
    takes the cut and puts all vertices in R on the sink side (removing them from cut) and places the vertices in C on the cut side (adding them to cut)
    '''
    cut = list(cut)
    for v in R:
        if v in cut:
            cut.remove(v)
    for v in C:
        cut.append(v)
    return tuple(cut)


def check_pfe_intersection(H,pfe,R):
    '''
    Finds if R intersects the tail of pfe but not all of the tail. Returns the tail vertices of pfe that R does not intersect and a boolean marking if R
    reaches some but not all of pfe's tail
    '''
    doesintersect = False
    tail = H.get_hyperedge_tail(pfe)
    reachedtailset = []
    for v in R:
        if v in tail:
            doesintersect = True
            reachedtailset.append(v)

    alltailreached = True
    for v in tail:
        if v not in reachedtailset:
            alltailreached = False
            break
    if alltailreached == True:
        doesintersect = False
    return reachedtailset,doesintersect


def find_reachable_vertices(H,edgesubset,complement):
    '''
    Finds the vertices in H that are reachable from the edges in edgesubset but not from the edges in complement
    Could make a new vertex in the hypergraph that is a new source vertex that just has one edge from it to all the vertices in the heads
    of all the edges in edgesubset, then call b visit on that new vertex and then delete the new vertex.
    '''
    newsource = 'fakesource'
    H.add_node(newsource)
    l = set()
    for e in edgesubset:
        for v in H.get_hyperedge_head(e):
            l.add(v)
        #for v in H.get_hyperedge_tail(e):
            #l.add(v)
    H.add_hyperedge([newsource],list(l))
    R = b_visit(H,newsource)

    newsource2 = 'fakecomplement'
    H.add_node(newsource2)
    l2 = set()
    for e in complement:
        for v in H.get_hyperedge_head(e):
            l2.add(v)
        #for v in H.get_hyperedge_tail(e):
            #l2.add(v)
    H.add_hyperedge([newsource2],list(l2))
    Rcomplement = b_visit(H,newsource2)

    H.remove_node(newsource)
    H.remove_node(newsource2)

    R = R[0]
    Rcomplement = Rcomplement[0]
    R.remove(newsource)
    Rcomplement.remove(newsource2)

    for v in Rcomplement:
        if v in R:
            R.remove(v)
    
    return R


def enumerate_edge_subsets(E):
    '''
    Enumerates the edge subsets of E up to a certain cardinality. Returns these as tuples of (edge subset, complement)
    '''
    E = set(E)
    boldF = []
    #TODO: change 
    MAX_CARD = 3 #maybe try not limiting the cardinality
    for i in range(1,MAX_CARD):
        sets = itertools.combinations(E,i)
        for s in sets:
            s = set(s)
            complement = E - s
            boldF.append((s,complement))
    return boldF

def find_crossing_edges(H,cj,normalcuts):
    '''
    Finds the edges that cross the cut cj
    '''
    E = []
    #first make sure the entire tail is in cj
    Cj = normalcuts[cj][0]
    for e in H.hyperedge_id_iterator():
        tailin = True
        for v in H.get_hyperedge_tail(e):
            if v not in Cj:
                tailin = False
                break
        #Then make sure that at least one vertex in the head is not in cj
        if tailin == True:
            for v in H.get_hyperedge_head(e):
                if v not in Cj:
                    E.append(e)
                    break

    return E

def find_head_cuts(H,taildistancelist,node_dict):
    '''
    Finds the normal cuts which just union the tails of edges with the same or smaller tail distance for all tail distance values in the sorted tail distance list.
    '''
    C = []
    V = set(['SUPERSOURCE'])
    distances = [taildistancelist[a] + 1 for a in taildistancelist]
    distances.append(0)
    distances = sorted(set(distances))
    distances.append(1000000000000)
    #print('distances',distances)
    vertexdistancelist = {}
    for v in H.get_node_set():
        mindist = 100000
        for e in H.get_backward_star(v):
            if taildistancelist[e] + H.get_hyperedge_weight(e) < mindist:
                mindist = taildistancelist[e] + H.get_hyperedge_weight(e)
        vertexdistancelist[v] = mindist
    vertexdistancelist['SUPERSOURCE'] = 0
            
    #print(vertexdistancelist)
    for d in distances:
        if vertexdistancelist['SUPERTARGET'] < d:
            break
        for v in H.get_node_set():
            if v in vertexdistancelist and vertexdistancelist[v] < d and v != 'SUPERTARGET':
                V.add(v)
        #print([node_dict[v] for v in V],d)
        C.append(tuple(V))

    return C

def find_normal_cuts(H,taildistancelist,node_dict):
    #TODO ultimately we need to prune the graph before we will be able to find the normal cuts that don't include anything that is not backward recoverable and reachable.
    '''
    Finds the normal cuts which just union the tails of edges with the same or smaller tail distance for all tail distance values in the sorted tail distance list.
    '''
    C = []
    V = set(['SUPERSOURCE'])
    distances = [taildistancelist[a] for a in taildistancelist]
    distances = sorted(set(distances))
    distances.append(1000000000000)
    #print('distances',distances)
    for d in distances:
        for e in H.hyperedge_id_iterator():
            if e in taildistancelist and taildistancelist[e] < d:
                for v in H.get_hyperedge_tail(e):
                    V.add(v)
        #print([node_dict[v] for v in V],d)
        C.append((tuple(V),d))

    return C


def find_pseudofinal_edges(H):
    '''
    Finds the pseudofinal edges, which are edges with multiple vertices in their tail set
    '''
    pseudofinaledges = []
    for e in H.hyperedge_id_iterator():
        if len(H.get_hyperedge_tail(e)) > 1:
            pseudofinaledges.append(e)
    return pseudofinaledges

def find_smallest_inedge(H,pfe,taildistancelist):
    '''
    Finds the inedge to the pseudofinal edge pfe with smallest d(Tail(e)) + weight(e)
    '''
    inedges = set()
    for v in H.get_hyperedge_tail(pfe):
        for f in H.hyperedge_id_iterator():
            if v in H.get_hyperedge_head(f):
                inedges.add(f)

    inedges = list(inedges)    
    e = inedges[0]
    edist = taildistancelist[e] + weight(H,[e])
    for f in inedges:
        if taildistancelist[f] + weight(H,[f]) < edist:
            e = f
            edist = taildistancelist[f] + weight(H,[f])
    return e

def find_cutgap_index(cutgapval,normalcuts):
    '''
    Finds the indices of the cutgap cuts that we will need to go back and augment
    cutgapval is a tuple the tail distance of the first cut and the last cut, so the start is the first normal cut that has that as its tail distance
    '''
    start = 0
    for i in range(len(normalcuts)):
        if cutgapval[0] <= normalcuts[i][1]:
            start = i
            break
    end = 0
    for i in range(len(normalcuts)):
        if cutgapval[1] <= normalcuts[i][1]:
            end = i
            break
    return start,end


def convert_vertex_cuts_to_edge_cuts(H,cuts,node_dict):
    edgecuts = []
    for cut in cuts:
        #print('converting cut',cut)
        if isinstance(cut,str):
            print(cut)
        if 'SUPERTARGET' in cut:
            cut.remove('SUPERTARGET')
            if 'SUPERSOURCE' not in cut:
                print('cut had supertarget and not supersource')
                cut.append('SUPERSOURCE')
            else:
                print('cut had supertarget')
        elif 'SUPERSOURCE' not in cut:
            print('cut had not supersource',len(cut))
            cut.append('SUPERSOURCE')
        edgecut = []
        for edge in H.hyperedge_id_iterator():
            tailin = True
            for v in H.get_hyperedge_tail(edge):
                if v not in cut:
                    tailin = False
                    break
            if tailin == True:
                for v in H.get_hyperedge_head(edge):
                    if v not in cut:
                        edgecut.append(edge)
                        break
                        
        if len(edgecut) == 0:
            print('no edges crossing cut')
            #print(('this cut did not have any edges crossing it',[node_dict[v] for v in cut]))
        #print('converted into', edgecut)
        edgecuts.append(tuple(edgecut))
    #print(len(edgecuts),'lenedgecuts')
    return set(edgecuts)



#######################################################################################################################################################################

##########################    New CUT STUFF 10/21 #################



def augment_head_cuts(H,headcuts,tailpathlist):
    '''
    returns the augmented head cuts (will be in the correct order, but internally they will be stored as a list of lists to allow for multiple versions of the same head cut.)
    '''
    print('headcutlens')
    for c in headcuts:
        print(len(c))
    #start by converting the headcuts into a list of lists
    headcuts = [[c] for c in headcuts]
    #print('loop headcuts',headcuts)
    extracuts = []
    #start by getting a list of the head cuts that each edge crosses?
    edges_crossing_headcuts, headcuts_crossed_by_edges,non_crossing_edges = find_cut_crossings(H,headcuts)
    #Loop through the head cuts looking for edges that crossed previous head cuts
    for i in range(len(headcuts)):
        headcut = headcuts[i][0]
        multiedges_crossing_i = []
        for e in headcuts_crossed_by_edges[i]:
            if len(edges_crossing_headcuts[e]) > 1 and i != min(edges_crossing_headcuts[e]):
                multiedges_crossing_i.append(e)
        #Call find_top_cut() on each of the edges that are multi-crossing
        top_cuts = []
        for e in multiedges_crossing_i:
            top_cut,f = find_top_cut(H,e,tailpathlist[e],headcuts,i,non_crossing_edges)
            if top_cut != None:
                top_cuts.append((top_cut,e,f))
        top_cuts.sort()
        #call process_multi_crossing_edge() on each multi-crossing edge, in order of their top cut

        for top_cut,e,f in top_cuts:
            topcuts,bottomcut = process_multi_crossing_edge(H,e,tailpathlist[e],headcuts[top_cut][0],headcut,f,non_crossing_edges)
            if topcuts != None:
                topcuts2 = None
                bottomcut2 = None
                #print('before',headcuts[top_cut][0], headcut)
                #print('after ',topcuts[0], bottomcut)
                #print(topcuts)
                extracuts = extracuts + list(headcuts[top_cut])
                extracuts.append(headcut)
                #extracuts.append(headcuts[top_cut])
                headcuts[i][0] = bottomcut
                if len(headcuts[top_cut]) > 1:
                    topcuts2,bottomcut2 = process_multi_crossing_edge(H,e,tailpathlist[e],headcuts[top_cut][1],headcut,f,non_crossing_edges)
                if len(topcuts) > 2:
                    #pick one that e crosses, which can just be returned as the first one
                    if topcuts2 != None:
                        headcuts[top_cut] = [topcuts[1],topcuts2[0]]
                        extracuts = extracuts + topcuts[2:] + topcuts2[1:] + [topcuts[0]]
                    else:
                        headcuts[top_cut] = topcuts[0:2]
                        extracuts = extracuts + topcuts[2:]
                else:
                    #print('replaced the old topcuts')
                    if topcuts2 != None:
                        headcuts[top_cut] = [topcuts[0],topcuts2[0]]
                    else:
                        headcuts[top_cut] = topcuts
                    #print('non-crossing',non_crossing_edges)
                edges_crossing_headcuts, headcuts_crossed_by_edges,non_crossing_edges = find_cut_crossings(H,headcuts)
                #print('extracuts',extracuts)
            
    headcuts = headcuts + extracuts
    return headcuts

def find_cut_crossings(H,headcuts):
    edges_crossing_headcuts = {}
    headcuts_crossed_by_edges = {i:[] for i in range(len(headcuts))}
    for i,headcut in enumerate(headcuts):
        for e in H.hyperedge_id_iterator():
            tail_in_source = True
            head_in_sink = False
            for v in H.get_hyperedge_tail(e):
                if v not in headcut[0]:
                    tail_in_source = False
                    break
            for v in H.get_hyperedge_head(e):
                if v not in headcut[0]:
                    head_in_sink = True
                    break
            if tail_in_source == True and head_in_sink == True:
                if e in edges_crossing_headcuts:
                    edges_crossing_headcuts[e].append(i)
                else:
                    edges_crossing_headcuts[e] = [i]
                headcuts_crossed_by_edges[i].append(e)
    
    non_crossing_edges = []
    #print('headcuts crossing',headcuts_crossed_by_edges)
    #print('edges crossing',edges_crossing_headcuts)
    for e in H.hyperedge_id_iterator():
        if e not in edges_crossing_headcuts:
            non_crossing_edges.append(e)
    return edges_crossing_headcuts, headcuts_crossed_by_edges,non_crossing_edges


def find_top_cut(H,e,P,headcuts,e_head_cut_index,non_crossing_edges):
    '''
    Takes in a multi-crossing edge e, and its path to its tail P and the head cuts
    returns the top cut for e, along with the edge f that is in P and is the non-crossing edge or an extra crossing
    '''
    #TODO Figure out why it sometimes cannot find a top cut. Is it likely due to non-crossing edges
    #TODO Maybe start by knowing which edges in the hyperpath to e are non-crossing, then we can ensure that they start completely only the source side of the bottom cut. Then we know we have reached the top cut when any part of their tail is moved to the sink side (because that shows that they will be crossed traversing up the path before you reach an edge actually crossing the top cut?)
    #print('finding top cut for ', headcuts[e_head_cut_index][0],H.get_hyperedge_tail(e),H.get_hyperedge_head(e))
    #print('non-crossing edges',non_crossing_edges)
    P_processed = []
    #print('path is')
    #for f in P:
        #print(f,H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
    topcut = headcuts[e_head_cut_index - 1]
    top_edge = None
    #loop backward through the head cuts while backward traversing the path
    for i in range(e_head_cut_index,-1,-1):
        curr_headcut = headcuts[i][0]
        #print(curr_headcut,'curr')
        P_edges_crossing = []
        P_noncrossing = []
        #print(P_processed)
        for f in P:
            if f not in P_processed or f in non_crossing_edges:
                f_tail_in_source = True
                f_head_in_sink = False
                for v in H.get_hyperedge_tail(f):
                    if v not in curr_headcut:
                        f_tail_in_source = False
                for v in H.get_hyperedge_head(f):
                    if v not in curr_headcut:
                        f_head_in_sink = True
                        if f_tail_in_source == True:
                            P_edges_crossing.append(f)
                            P_processed.append(f)
                            break
                if f_head_in_sink == True or f_tail_in_source == False:
                    #print(f,'bottom')
                    #If its head is on the sink side, then for the bottom cut, its tail is on the source side, so pulling over its tail vertex would make it start crossing the cut before we got to an edge that crosses this cut.
                    if f in non_crossing_edges:
                        P_noncrossing.append(f)
                        #print('non-crossing edge',H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
        if len(P_noncrossing) > 0:
            #print('non-crossing',P_noncrossing[0])
            topcut = i
            top_edge = P_noncrossing[0]
            return topcut, top_edge
        elif len(P_edges_crossing) > 1:
            #print('extra-crossing',P_edges_crossing)
            topcut = i
            top_edge = P_edges_crossing[0]
            return topcut, top_edge
        #print(curr_headcut,P_edges_crossing)
    print('error, did not find top cut for {} {} {} with tail len {}'.format(e,P,e_head_cut_index,len(H.get_hyperedge_tail(e))))
    print('non-crossing P edges',[f for f in non_crossing_edges if f in P])
    for i in range(e_head_cut_index):
        curr_headcut = headcuts[i][0]
        noncrossp = [f for f in non_crossing_edges if f in P]
        print(noncrossp,'noncrossp',i,len(curr_headcut))
        for f in noncrossp:
            f_tail_in_source = True
            f_head_in_sink = False
            for v in H.get_hyperedge_tail(f):
                if v not in curr_headcut:
                    f_tail_in_source = False
            for v in H.get_hyperedge_head(f):
                if v not in curr_headcut:
                    f_head_in_sink = True
            if f_head_in_sink == True or f_tail_in_source == False:
                print(f,i)

                    
            
    #return None, None
    return e_head_cut_index - 1, non_crossing_edges[0]
                

    #stop when either the procedure passes a non-crossing edge, or reaches a cut with multiple crossings
    #return the non-crossing edge or the extra crossing


def process_multi_crossing_edge(H,e,P,topcut,bottomcut,f,non_crossing_edges):
    '''
    returns the list of head cuts that have been augmented for e
    takes in the edge that is being processed, a path to the tailset of that edge, the headcuts, its top cut, its bottom cut, the edge you are making cross the bottom cut
    '''

    #start by finding the quasi-reachable in-edges to e from the top edge f, called e_in
    #Find the edges in P that are quasi-reachable from f by crawling down the path
    vertices_reached = set(H.get_hyperedge_head(f))
    #print('vertices_reached',vertices_reached)
    e_tail_set = set(H.get_hyperedge_tail(e))
    #print('tailset e',e_tail_set)
    reached_edges = [f]
    e_in = f
    #find the in-edges to e that are reachable from f (one is picked arbitrarily when multiple exist)
    while len(e_tail_set.intersection(vertices_reached)) == 0:
        if len(reached_edges) == 0:
            return None, None
        e_in = reached_edges[0]
        for v in H.get_hyperedge_head(e_in):
            vertices_reached.add(v)
            for g in H.get_forward_star(v):
                if g in P and g not in reached_edges:
                    reached_edges.append(g)
        reached_edges.remove(e_in)
    v_in = 0
    #print(len( set(H.get_hyperedge_head(e_in)).intersection(set(H.get_hyperedge_tail(e)))))
    for v in set(H.get_hyperedge_head(e_in)).intersection(set(H.get_hyperedge_tail(e))):
        v_in = v
        v_count = 0
        for g in H.get_backward_star(v):
            if g in P:
                v_count += 1
        if v_count == 1:
            break
    if v_in == 0:
        print('error, v_in was not set')
    bottomcut = list(bottomcut)
    originalbottomcut = bottomcut
    bottomcut, E_top = upward_cut_augment(H,e,P,[],topcut,bottomcut,non_crossing_edges,v_in=v_in)
    E_bottom = []
    #for g in H.get_backward_star(v_in):
        #g_tail_in_source = True
        #for w in H.get_hyperedge_tail(g):
            #if w not in bottomcut:
                #g_tail_in_source = False
        #if g_tail_in_source == True:
            #E_bottom.append(g)

    #Then perform an upward pass from the vertex
    #bottomcut.remove(v_in)
    previous_E_top = set(E_top)
    previous_E_bottom = set(E_bottom)
    #print('beginning a thing')
    #while len(E_bottom) > 0:
    while len(E_top) > 0:
        E_bottom_to_remove = []
        E_top_to_remove = []
        #Need to keep track of what has already been in E_top and E_bottom so that we don't do the same things over and over again
        #bottomcut, E_top = upward_cut_augment(H,e,P,E_bottom,topcut,bottomcut,non_crossing_edges)
        topcut, E_bottom = downward_cut_augment(H,e,P,E_top,topcut,bottomcut,originalbottomcut)
        if 'SUPERSOURCE' not in bottomcut:
            print('supersource not in bottomcut')
        if 'SUPERTARGET' in topcut:
            print('supertarget in topcut')
        #print('E top',E_top,[(H.get_hyperedge_tail(h),H.get_hyperedge_head(h)) for h in E_top])
        for g in E_bottom:
            if g in previous_E_bottom:
                E_bottom_to_remove.append(g)
            else:
                previous_E_bottom.add(g)
        E_bottom = [b for b in E_bottom if b not in E_bottom_to_remove]
        #if len(E_top) > 0:
        if len(E_bottom) > 0:
            #topcut, E_bottom = downward_cut_augment(H,e,P,E_top,topcut,bottomcut,originalbottomcut)
            bottomcut, E_top = upward_cut_augment(H,e,P,E_bottom,topcut,bottomcut,non_crossing_edges)
            if 'SUPERSOURCE' not in bottomcut:
                print('supersource not in bottomcut')
            if 'SUPERTARGET' in topcut:
                print('supertarget in topcut')
            for g in E_top:
                if g in previous_E_top:
                    E_top_to_remove.append(g)
                else:
                    previous_E_top.add(g)
            E_top = [b for b in E_top if b not in E_top_to_remove]
        else:
            #E_bottom = []
            E_top = []
            break
    topcuts = [topcut]

    e_tail_in_headcut = True
    for v in H.get_hyperedge_tail(e):
        if v not in topcut:
            e_tail_in_headcut = False
            break

    if e_tail_in_headcut:
        #print('made alternatecuts \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n')
        #print(H.get_hyperedge_tail(e))
        #print(topcut)
        #print(topcuts)
        #make alternate versions
        for v in H.get_hyperedge_tail(e):
            alt_topcut = [w for w in topcut if w != v]
            #print('alt',alt_topcut)

            topcuts.append(alt_topcut)
        #print(topcuts)
    return topcuts,bottomcut

def upward_cut_augment(H,e,P,E_bottom,topcut,bottomcut,non_crossing_edges,v_in=None):
    '''
    returns the augmented bottom cut and a list of edges that cross the top cut that now need to be explored
    Takes in the edge that is being processed, the path to that edge, the set of edges that now need to be moved to the sink side of the cut, the top cut and the bottom cut.
    '''
    #print('starting upward cut augment')
    #print('E_bottom',E_bottom)
    #print('topcut',topcut)
    edges_to_explore = E_bottom
    E_top = set()
    edges_to_remove = []
    for f in edges_to_explore:
        #make sure that none of these edges cross the top cut
        #print('exploring {}'.format(f))
        f_tail_in_topcut = True
        for v in H.get_hyperedge_tail(f):
            if v not in topcut:
                f_tail_in_topcut = False
                break
        if f_tail_in_topcut == True:
            #print('{} had its tail in the topcut'.format(f))
            #This means that the edge in E_bottom already crosses E_top, so nothing more needs to be done
            E_top.add(f)
            edges_to_remove.append(f)
    edges_to_explore = [e for e in edges_to_explore if e not in edges_to_remove]

    bottomcut = list(bottomcut)

    if v_in != None:
        #we know the vertex that we are going to remove from the bottom cut, now we just need to figure out which edges that causes to cross
        bottomcut.remove(v_in)
        newedges = H.get_backward_star(v_in)

        #First find the edges that are now crossing
        for f in newedges:
            f_crosses_topcut = True
            for w in H.get_hyperedge_tail(f):
                if w not in topcut:
                    f_crosses_topcut = False
                    break
            if f_crosses_topcut == False:
                if f not in non_crossing_edges:
                    edges_to_explore.append(f)
            else:
                E_top.add(f)


    #print('E_bottom after trimming',edges_to_explore)
    #starting with the tailset of each edge in E_bottom, call find_smallest_new_edge_crossings() on each edge
    while len(edges_to_explore) > 0:
        #print('edges to explore',edges_to_explore)
        new_edges_to_explore = []
        for e in edges_to_explore:
            v,v_count,newedges = find_smallest_new_edge_crossings(H,e,bottomcut)
            #print('edge {} with tail {} and head {} removes {} because it has {} new edges {} crossing the bottom cut {}'.format(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e),v,v_count,newedges,bottomcut))
            if v != None:
                bottomcut.remove(v)
                for f in newedges:
                    f_crosses_topcut = True
                    for w in H.get_hyperedge_tail(f):
                        if w not in topcut:
                            f_crosses_topcut = False
                            break
                    if f_crosses_topcut == False:
                        if f not in non_crossing_edges:
                            new_edges_to_explore.append(f)
                    else:
                        E_top.add(f)
        edges_to_explore = new_edges_to_explore
                    
    return bottomcut, list(E_top)
            

    #Then remove those edges from the bottom cut and see what new edges cross the bottom cut and repeat checking if the new edges crossing are crossing the top cut

def find_smallest_new_edge_crossings(H,e,bottomcut):
    '''
    returns the vertex to be taken out of the sink of the cut that makes the fewest new edges cross the cut, along with a list of the new edges that will cross the cut

    '''

    #loop through the vertices v in the tail of e, seeing which edges with v in their head have their tail in bottomcut
    best_vertex = None
    best_vertex_count = 999999
    new_edges = []
    for v in H.get_hyperedge_tail(e):
        if v in bottomcut:
            v_backstar = H.get_backward_star(v)
            v_new_edges = []
            for f in v_backstar:
                tail_in_source = True
                for w in H.get_hyperedge_tail(f):
                    if w not in bottomcut:
                        tail_in_source = False
                if tail_in_source == True:
                    v_new_edges.append(f)
            v_count = len(v_new_edges)
            if v_count < best_vertex_count:
                best_vertex = v
                best_vertex_count = v_count
                new_edges = v_new_edges
    return best_vertex,best_vertex_count,new_edges


#def downward_cut_augment(H,e,P,E_top,topcut,bottomcut):
def downward_cut_augment(H,e,P,E_top,topcut,bottomcut,originalbottomcut):
    '''
    returns the augmented top cut and a list of edges that cross the bottom cut that now need to be removed from the bottom cut
    Takes in the edge that is being processed, the path to that edge, the set of edges that now need to be moved to the source side of the cut, the top cut and the bottom cut.
    '''
    #print('starting downward cut augment')
    #print('E_top')
    #for h in E_top:
        #print(h,H.get_hyperedge_tail(h),H.get_hyperedge_head(h))
    E_bottom = set()
    edges_to_explore = E_top
    
    explore_edges_to_remove = []
    for g in edges_to_explore:
        tail_in_source = True
        head_crosses_bottom = False
        for v in H.get_hyperedge_tail(g):
            if v not in originalbottomcut:
                tail_in_source = False
                break
        if tail_in_source == True:
            for w in H.get_hyperedge_head(g):
                if w not in originalbottomcut:
                    head_crosses_bottom = True
                    explore_edges_to_remove.append(g)
                    break
    edges_to_explore = [g for g in edges_to_explore if g not in explore_edges_to_remove]
        
    topcut = list(topcut)
    #Until the edges that now cross the cut previously crossed the bottom cut or there aren't any more
    while len(edges_to_explore) > 0:
        new_edges_to_explore = []
        for g in edges_to_explore:
            forward_stars = []
            for v in H.get_hyperedge_head(g):
                if v not in topcut:
                    topcut.append(v)
                    for f in H.get_forward_star(v):
                        if f not in forward_stars:
                            forward_stars.append(f)
            for f in forward_stars:
                if f != e:
                    tail_in_source = True
                    tail_in_bottom_source = True
                    head_crosses_bottom = False
                    tail_in_og_bottom_source = True
                    head_crosses_og_bottom = False
                    head_crosses_top = False
                    for w in H.get_hyperedge_tail(f):
                        if w not in topcut:
                            tail_in_source = False
                        if w not in bottomcut:
                            tail_in_bottom_source = False
                        if w not in originalbottomcut:
                            tail_in_og_bottom_source = False

                    for w in H.get_hyperedge_head(f):
                        if w not in bottomcut:
                            head_crosses_bottom = True 
                        if w not in originalbottomcut:
                            head_crosses_og_bottom = True 
                        if w not in topcut:
                            head_crosses_top = True
                    if (tail_in_bottom_source == True and head_crosses_bottom == True) or (tail_in_og_bottom_source == True and head_crosses_og_bottom == True):
                        #if it crosses the bottomcut
                        E_bottom.add(f)     
                    elif tail_in_source == True and head_crosses_top == True:
                        #if it crosses the top cut but not the bottom cut
                        new_edges_to_explore.append(f)
        edges_to_explore = new_edges_to_explore
    return list(topcut),list(E_bottom)

    
    #Put their head sets onto the source side of the cut

#######################################################################################################################################################################

#######################################################################################################################################################################
'''
test case:
vertices will be added automatically, so we just need to define the edges as tail and head sets
1 to 2
2 to 3
2 to 4
3 to 5
4 to 6
5,6 to 7

test case 2:
We want multiple paths in this case, with multiple edges with different top cuts but that interact
Also want a path to the entire tail set of e
'''
#######################################################################################################################################################################

if __name__ == '__main__':
    H = DirectedHypergraph()
    hyperedges = [
                    (['SUPERSOURCE'],[1]),
                    ([1],[2]),
                    ([2],[3]),
                    ([2],[4]),
                    ([4],[5]),
                    ([4],[13]),
                    ([3],[6]),
                    ([3],[9]),
                    ([3],[14]),
                    ([3],[5,13]),
                    ([4],[6,14]),
                    ([8],['SUPERTARGET']),
                    ([7],[8]),
                    ([2],[9]),
                    ([9],[10]),
                    ([10],[11]),
                    ([11],[12]),
                    ([12],[15]),
                    ([15],[16]),
                    ([16],[17]),
                    ([17],[18]),
                    ([18],['SUPERTARGET']),
                    ([2],[19]),
                    ([19],[20]),
                    ([20],[21]),
                    #([21],[22]),
                    #([22],[23]),
                    #([23],[24]),
                    #([24],[25]),
                    #([25],[26]),
                    #([26],[5,6,13,14]),
                    ([21],[5,6,13,14]),
                    ([5,6,13,14],[7])]

    H.add_hyperedges(hyperedges)
    
    taildistancelist,H,returnpath,tailpathlist = tail_path_heuristic(H,'SUPERSOURCE','SUPERTARGET',return_paths=True)
    print([(H.get_hyperedge_tail(e),H.get_hyperedge_head(e)) for e in returnpath])
    headcuts = find_head_cuts(H,taildistancelist,{})
    augheadcuts = augment_head_cuts(H,headcuts,tailpathlist)
    print()
    print(headcuts)
    print()
    print(augheadcuts)
    print()
    flataugmentedcuts = []
    for c in augheadcuts:
        if isinstance(c[0],list) or isinstance(c[0],tuple):
            for b in c:
                flataugmentedcuts.append(b)
        else:
            flataugmentedcuts.append(c)
    edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H,headcuts,{})
    edgeaugmentedcuts = convert_vertex_cuts_to_edge_cuts(H,flataugmentedcuts,{})
    print(edgeaugmentedcuts)
    ROOTDIR = '/Users/skrieger/Documents/UofA/biopax-hypergraph/biopax-hypergraph/test'
    outprefix = '%s/results/reactome-my-acyclic' % (ROOTDIR)
    ilpname = '%s/results/reactome-my-heuristic.lp' % (ROOTDIR)
    make_shortest_cyclic_hyperpath_ilp(H,'SUPERSOURCE','SUPERTARGET',ilpname)
    run_LP_heuristics(H,H.get_node_set(),ilpname,outprefix,1,'SUPERTARGET','SUPERSOURCE',edgeheadcuts,edgeaugmentedcuts,targetname='SUPERTARGET')
