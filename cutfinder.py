import itertools
import time

from ilp import *
from pathheuristic import weight, tail_path_heuristic
from pathheuristic import *

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit,visit


#######################################################################################################################################################################

def find_cuts_run_heuristic(H,source,target,node_dict,return_dict,event):

    '''
    Runs the heuristic and then Finds the cuts based off the heuristic path's tail distance list.
    '''

    taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True,doublyreachablesubgraph=True,print_paths=True,event=event)
    if event.is_set():
        return

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


def convert_vertex_cuts_to_edge_cuts(H,cutlist):
    '''
    Takes cuts in a cutlist that are defined by vertices in H and converts them into cuts that are based on edges in H
    '''
    edgecuts = set()
    for cut in range(len(cutlist)):
        E = find_crossing_edges(H,cut,cutlist)
        edgecuts.add(tuple(E))
    return edgecuts
        

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


