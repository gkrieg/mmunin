from xml.dom import minidom
import sys
import cplex
import itertools
import heapq
import json
import time
from multiprocessing import Pool

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit,visit
from halp.utilities import directed_graph_transformations
from queue import Queue

EPSILON=0.0001
CONSTANT=1000000

NODE_DICT = {}

def tail_path_heuristic(H,source,target,node_dict={},best_in_edges=False,print_paths=False,return_paths=False,sinkrecover=False,timingstats=False,doublyreachablesubgraph=False):
    '''
    Finds the distance to the tail for each edge in the hypergraph H that is reachable from the source vertex and backwards-recoverable from the target
    Returns a list of the tail distance for each edge
    ''' 
    begin = time.time()
    reachableedges,edgedict,taildistancelist,heap,reachableedgecounter,reachedtable,entry_finder,counter,H = initialize(H,source,target,node_dict=node_dict,alreadydoublyreachable=doublyreachablesubgraph)
    afterinitialization = time.time()
    print(('initialization took {0}'.format(afterinitialization - begin)))
    print('done initializing, starting while loop')
    tailpathslist = {}

    returnpath = []
    edgeprocesstime = 0
    while reachableedgecounter > 0:
        if reachableedgecounter % 200 == 0:
            print('number of edges still to process: {}'.format(reachableedgecounter),flush=True)
        e,epriority = pop_node(heap,entry_finder)
        
        if e in reachableedges and edgedict[e]['isremoved'] == False:
            reachableedgecounter -= 1
            #reachableedges.remove(e) #possibly don't need this
            edgedict[e]['isremoved'] = True
            P,inedges = recover_short_hyperpath(H,e,best_in_edges,edgedict,source,taildistancelist)

            if target in H.get_hyperedge_head(e):
                returnpath = print_out_sthyperpath(P,returnpath,print_paths,H,node_dict)

            taildistancelist[e] = weight(H,P) - weight(H,[e])
            tailpathslist[e] = P
            edgedict[e]['bestinedges'] = inedges
            edgestoprocess = findedgestoprocess(e,H,reachedtable,edgedict)
            #print('edges to process',edgestoprocess)
            edgeprocesstimebefore = time.time()
            '''
            for f in edgestoprocess:
                edgedict[f]['candidateinedges'].append(e)
                if f in entry_finder and edgedict[f]['isremoved'] == False:
                    fpath,_ = recover_short_hyperpath(H,f,best_in_edges,edgedict,source,taildistancelist)
                    edgepriority = weight(H,fpath)
                    add_node(heap,f,edgepriority,counter,entry_finder)
                elif edgedict[f]['tailcount'] == 0 and f not in entry_finder:
                    fpath,_ = recover_short_hyperpath(H,f,best_in_edges,edgedict,source,taildistancelist)
                    edgepriority = weight(H,fpath)
                    add_node(heap,f,edgepriority,counter,entry_finder)
            '''
            #NOTE: new parallel section
            edgestoupdate = []
            for f in edgestoprocess:
                edgedict[f]['candidateinedges'].append(e)
                if (f in entry_finder and edgedict[f]['isremoved'] == False) or (edgedict[f]['tailcount'] == 0 and f not in entry_finder):
                    edgestoupdate.append(f)
            with Pool(32) as pool:
                priorities = pool.starmap(process_edge,zip(edgestoupdate, itertools.repeat(H), itertools.repeat(best_in_edges), itertools.repeat(edgedict), itertools.repeat(source), itertools.repeat(taildistancelist)))
            for i,f in enumerate(edgestoupdate):
                edgepriority = priorities[i]
                add_node(heap,f,edgepriority,counter,entry_finder)
            #NOTE: end new parallel section

            if timingstats == True:
                print('time to process edges',edgeprocesstimeafter - edgeprocesstimebefore)
            edgeprocesstimeafter = time.time()
            edgeprocesstime += edgeprocesstimeafter - edgeprocesstimebefore
        else:
            print(('edge {0} was not in the reachable set but still in the hypergraph'.format(e)))

    #print(taildistancelist)
    if sinkrecover == True:
        recoveredlist = set(H.get_backward_star(target))
        recoveredlistall = set(H.get_backward_star(target))
        print(len(H.get_backward_star(target)),'lelelel')
        for f in H.get_backward_star(target):
            a,b,c = recover(H,f,'full',edgedict)
            a2,b2,c2 = recover(H,f,'all',edgedict)
            fpath,_ = trim(a2,b2,c2,source,taildistancelist)

            print(len(b2),'lenb2')
            print('pathlength',weight(H,fpath))
            for bb in b:
                recoveredlist.add(bb)
            for bb2 in b2:
                recoveredlistall.add(bb2)
        
    end = time.time()
    print(('while loop took {0}'.format(end - afterinitialization)))
    print(('whole path heuristic took {0}'.format(end - begin)))
    print('time spent processing edges {}'.format(edgeprocesstime))
    if return_paths == True:
        return taildistancelist,H,returnpath,tailpathslist
    elif sinkrecover == True:
        return taildistancelist,H,returnpath,recoveredlist,recoveredlistall
    else:
        return taildistancelist,H,returnpath

def process_edge(f,H,best_in_edges,edgedict,source,taildistancelist):
    fpath,_ = recover_short_hyperpath(H,f,best_in_edges,edgedict,source,taildistancelist)
    edgepriority = weight(H,fpath)
    return edgepriority

def print_out_sthyperpath(P,returnpath,print_paths,H,node_dict):
    print('Path from Trim')
    if len(P) < len(returnpath) or returnpath == []:
        returnpath = P
    if print_paths == True:
        for p in P:
            tail = []
            unreadabletail = []
            for v in H.get_hyperedge_tail(p):
                if len(node_dict) > 0:
                    tail.append(node_dict[v])
                    unreadabletail.append(v)
                else:
                    tail.append(H.get_node_attribute(v,'label'))
                    unreadabletail.append(v)
            head = []
            unreadablehead = []
            for v in H.get_hyperedge_head(p):
                if len(node_dict) > 0:
                    head.append(node_dict[v])
                    unreadablehead.append(v)
                else:
                    head.append(H.get_node_attribute(v,'label'))
                    unreadablehead.append(v)

            print((p,tail,head))
            #print(('unreadable',p,unreadabletail,unreadablehead))
    print(('with weight weight{0}\n'.format(weight(H,P))))
    return returnpath

def recover_short_hyperpath(H,f,which_inedges,edgedict,source,taildistancelist,timingstats=False):

    if timingstats == True:
        recoverbefore = time.time()

    if which_inedges == False:
        a,b,c = recover(H,f,'full',edgedict)
    else:
        a,b,c = recover(H,f,'short',edgedict)

    if timingstats == True:
        recoverafter = time.time()
        print('recover time',recoverafter - recoverbefore)

    fpath,inedges = trim(a,b,c,source,taildistancelist)

    if timingstats == True:
        trimafter = time.time()
        print('trim time',trimafter - recoverafter)

    return fpath,inedges


def heuristic_path(H,source,target):
    '''
    Does some heuristic method to find the shortest path so that we can use that path to seed the ILP with cuts
    This one in particular uses a modified version of Dijkstra's algorithm, which will approximately calculate the shortest path.

    '''
    entry_finder = {}               # mapping of tasks to entries
    counter = itertools.count()     # unique sequence count
    priorityq = []
    add_node(priorityq,source,0,counter,entry_finder)
    processed_nodes = set()
    node_dict = {} #has a tuple of (priority,set of edges used to get there)
    node_dict[source] = (0,[])
    curnode = source
    while not emptypq(priorityq,entry_finder) and curnode != target:
        curnode,curpriority = pop_node(priorityq,entry_finder)
        #print curnode,curpriority
        for edge in H.hyperedge_id_iterator():
            priority = 0
            edgelist = []
            for node in H.get_hyperedge_tail(edge):
                if node not in processed_nodes and node != curnode:
                    break
                priority += node_dict[node][0] #This is the heuristic weight where we sum the distance to all vertices in the tail
                edgelist += node_dict[node][1]
            #If we get to this point, the edge can be traversed
            else:
                for node in H.get_hyperedge_head(edge):
                    if node not in processed_nodes and (node not in node_dict or priority + H.get_hyperedge_weight(edge) < node_dict[node][0]):
                        add_node(priorityq,node,priority + H.get_hyperedge_weight(edge),counter,entry_finder)
                        node_dict[node] = (H.get_hyperedge_weight(edge) + priority,edgelist + [edge])

        processed_nodes.add(curnode)

    return node_dict[curnode][1]



def findedgestoprocess(e,H,reachedtable,edgedict):
    '''
    Finds the edges whose tails intersect with e's head. Also sets e to their inedge lists and decrements their reached counter in the reachedtable
    '''
    F = set()
    #print('finding edges for {} which has tail {} and head {}'.format(e,[NODE_DICT[v.strip()] for v in H.get_hyperedge_tail(e) if v in NODE_DICT],[NODE_DICT[v] for v in H.get_hyperedge_head(e)]))
    for v in H.get_hyperedge_head(e):
        #print(v)
        for f in H.get_forward_star(v):
            #edgedict[f]['candidateinedges'].append(e)
            #print('found edge {} with {} it its tail'.format(f,v))
            if v not in reachedtable[f]:
                reachedtable[f].append(v)
                #if f == 'e11':
                    #print(edgedict[f]['tailcount'],e)
                edgedict[f]['tailcount'] -= 1
            if v not in reachedtable[f] or edgedict[f]['isremoved'] == False:
                F.add(f)
                #print('adding {} to F for {}'.format(f,e))
            #else:
                #print('{} was already in the reached table for {}'.format(v,f))

    return list(F)


def weight(H,F):
    '''
    Takes as input a set of edges and returns the sum of the weight of the edges
    '''
    w = 0
    for f in F:
        w += H.get_hyperedge_weight(f)
    return w

def recover(H,e,flag,edgedict):
    '''
    finds the set of edges that are recovered recursively for each edge starting with the inedge list of e. Returns this set of edges
    The flag sets whether to use bestinedges or candidate inedges for every edge except e, or to use all inedges for all edges
    '''
    Q = []
    Qindex = 0
    F = set()
    #if e == 'e11' and flag == 'full':
        #print('printing es candidates')
        #print(edgedict[e]['candidateinedges'])
    if flag == 'all':
        #need to find the set of all inedges for all edges
        #newedgedict = {}
        #for f in H.hyperedge_id_iterator():
            #inedges = set()
            #for v in H.get_hyperedge_tail(f):
                #edgelist = H.get_backward_star(v)
                #for g in edgelist:
                    #inedges.add(g)
                ##for g in H.hyperedge_id_iterator():
                    ##if v in H.get_hyperedge_head(g):
                        ##inedges.add(g)
            ##inedges = H.get_backward_star(f)
            #newedgedict[f] = list(inedges)
        #for f in newedgedict[e]:
        for v in H.get_hyperedge_tail(e):
            for f in H.get_backward_star(v):
                if f not in F:
                    Q.append(f)
                    F.add(f)
    if flag != 'all':
        for f in edgedict[e]['candidateinedges']:
            Q.append(f)
            F.add(f)
    while Qindex < len(Q):
        #process the next edge in the Q
        f = Q[Qindex]
        Qindex += 1
        if flag == 'all':
            for v in H.get_hyperedge_tail(f):
                for g in H.get_backward_star(v):
                    if g not in F:
                        F.add(g)
                        Q.append(g)
        elif flag == 'full':
            for g in edgedict[f]['candidateinedges']:
                if g not in F:
                    F.add(g)
                    Q.append(g)
        elif flag == 'short':
            for g in edgedict[f]['bestinedges']:
                if g not in F:
                    F.add(g)
                    Q.append(g)
        else:
            print('something wrong in recover with the flag')

    #if flag == 'all':
        #print('lenQ',Qindex)
    return H,F,e

def trim(H,F,e,source,taildistancelist,timingstats=False):
    '''
    Takes the edge list F, which is a super path from the source vertex to e and repeatedly removes an edge and checks reachability for each edge in F
    '''
    #should first see if there is reachability before pruning any edge, if so, it returns no solution
    #start by sorting the edges into the order you are going to consider them
    if timingstats == True:
        starttime = time.time()
    if e[0] != 'e':
        reachingset = [e]
    else:
        reachingset = H.get_hyperedge_tail(e)
        
    sortededges = []
    #tdistlist = [taildistancelist[t] for t in taildistancelist if taildistancelist[t] != 'inf' and taildistancelist[t] > 0]
    #if len(tdistlist) > 0:
        #print(tdistlist)
    for f in F:
        sortededges.append((taildistancelist[f],f))
    sortededges.sort(reverse=True)
    justedges = [s[1] for s in sortededges]
    if timingstats == True:
        sorttime = time.time()
        print('time to sort',sorttime - starttime)
    H2 = DirectedHypergraph()
    #H2.add_nodes(H.get_node_set())
    H2.add_node(source)
    H2edgeids = {}
    Hedgeids = {}
    for f in justedges:
        newf = H2.add_hyperedge(H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
        H2edgeids[f] = newf
        Hedgeids[newf] = f
    if timingstats == True:
        buildtime = time.time()
        print('time to build alternate hypergraph',buildtime-sorttime)
    if check_reachability(H2,reachingset,source) == False:
        print('invalid edge set given to trim. Sink not reachable from source')
        print('number of edges in new hypergraph',len(H2.get_hyperedge_id_set()))

    #nodeset,_,__,___ = b_visit(H2,source)
    if timingstats == True:
        bvisittime = time.time()
        print('time for b visit',bvisittime-buildtime)
    #isereached = True
    #for v in reachingset:
        #if v not in nodeset:
            #isereached = False
    #if isereached == False:
        #print('invalid edge set given to trim. Sink not reachable from source')
        #print(nodeset)
    F2 = []
    tailedges = []
    forbiddenset = set()
    for f in justedges:
        forbiddenset.add(H2edgeids[f])
        if check_reachability(H2,reachingset,source,forbiddenedges=forbiddenset) == False:
            forbiddenset.remove(H2edgeids[f])
            F2.append(f)
            #if e[0] == 'e':
                #for v in H.get_hyperedge_tail(e):
                    #if v in H.get_hyperedge_head(f):
                        #tailedges.append(f)

    tailedges = set()
    if e[0] == 'e':
        F2.append(e)
        #collecting the edges going into the tail
        for v in H.get_hyperedge_tail(e):
            for f in H2.get_backward_star(v):
                if f not in forbiddenset:
                    tailedges.add(Hedgeids[f])
        #nodeset,_,__,___ = b_visit(H2,source)
        #isereached = True
        #for v in reachingset:
            #if v not in nodeset:
                #isereached = False
        #if isereached == False:
            #This means we cannot remove f
            #newf = H2.add_hyperedge(H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
            #H2edgeids[f] = newf
            #F2.append(f)
            #if e[0] == 'e':
                #for v in H.get_hyperedge_tail(e):
                    #if v in H.get_hyperedge_head(f):
                        #tailedges.append(f)
    if timingstats == True:
        whiletime = time.time()
        print('time in while loop', whiletime - bvisittime)
    return F2,list(tailedges)

def check_reachability(H,reachset,source,forbiddenedges=set()):
    
    '''
    Basically B_visit but stops and returns True when all the vertices are reached

    '''
    reachsetreachedcounter = len(reachset)
    if source in reachset:
        reachsetreachedcounter -= 1
        if reachsetreachedcounter == 0:
            return True

    visited_nodes = set([source])
    k = {hyperedge_id: len(H.get_hyperedge_tail(hyperedge_id)) for hyperedge_id in H.get_hyperedge_id_set()}

    Q = Queue()
    Q.put(source)

    while not Q.empty():
        current_node = Q.get()
        # At current_node, we can traverse each hyperedge in its forward star
        for hyperedge_id in H.get_forward_star(current_node):
            # Since we're arrived at a new node, we increment
            # k[hyperedge_id] to indicate that we've reached 1 new
            # node in this hyperedge's tail
            if hyperedge_id not in forbiddenedges:
                k[hyperedge_id] -= 1
                # Traverse this hyperedge only when we have reached all the nodes
                # in its tail (i.e., when k[hyperedge_id] == |T(hyperedge_id)|)
                if k[hyperedge_id] == 0:
                    #Pe[hyperedge_id] = current_node
                    # Traversing the hyperedge yields the set of head nodes of
                    # the hyperedge; B-visit each head node
                    for head_node in H.get_hyperedge_head(hyperedge_id):
                        if head_node in visited_nodes:
                            continue
                        if head_node in reachset:
                            reachsetreachedcounter -= 1
                            if reachsetreachedcounter == 0:
                                return True
                        Q.put(head_node)
                        visited_nodes.add(head_node)
    return False

def get_doubly_reachable_graph(H,source,target,node_dict):
    reachableedges = findreachableandbackrecoverable(H,source,target)
    print(source,target)
    print('len reachable',len(reachableedges))
    numedges = 0
    H2 = DirectedHypergraph()
    if len(node_dict) == 0:
        for v in H.node_iterator():
            H2.add_node(v,H.get_node_attributes(v))
    for edge in reachableedges:
        if target in H.get_hyperedge_head(edge):
            print(edge,H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge))
        eid = H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge),weight =H.get_hyperedge_weight(edge))
        numedges += 1
    #H2.add_node('SUPERSOURCE',{'label': 'SUPERSOURCE'}) 
    for edge in H2.get_forward_star(source):
        forwardstarlist = []
        for v in H2.get_hyperedge_head(edge):
            if len(H2.get_forward_star(v)) > 0:
                forwardstarlist.append(v)
        H2.remove_hyperedge(edge)
        numedges -= 1
        if len(forwardstarlist) > 0:
            H2.add_hyperedge([source],forwardstarlist,weight=0)
            numedges += 1
        else:
            H2.add_hyperedge([source],[],weight=0)
            numedges += 1
    print('got doubly reachable graph. It has {} edges from {} reachableedges'.format(len(H2.get_hyperedge_id_set()),len(reachableedges)))
    print(numedges)
    return H2

def initialize(H,source,target,node_dict={},alreadydoublyreachable=False):
    '''
    Finds the set of reachable and backwards-recoverable edges, sets up the tail reached counters, the inedge lists,  and the heap pointer for each edge.
    nitializes the taillength for each edge, and the heap, and the reachable edge counter
    '''
    edgedict = {}

    #find reachable and backwards-recoverable edge list
    #node_dict = printremappededges(H)
    #print('found reachable edges')

    #trim H

    #print('building new hypergraph from doubly reachable set')
    if alreadydoublyreachable == False:
        reachableedges = findreachableandbackrecoverable(H,source,target)
        reachableedges.sort()
        print('initializing hypergraph with reachable edges. Len reachable edges is {}, len hyperedgeidset is {}'.format(len(reachableedges),len(H.get_hyperedge_id_set())))
        H2 = DirectedHypergraph()
        if len(node_dict) == 0:
            for v in H.node_iterator():
                H2.add_node(v,H.get_node_attributes(v))
        #H2.add_nodes(H.get_node_set())
        #print('done adding new nodes')
        #print('targetH',target, H.get_backward_star(target))
        #for backedge in H.get_backward_star(target):
            #if backedge in reachableedges:
                #print('backhere')
            #else:
                #print('notbackhere')
        for edge in reachableedges:
            #H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge),{'oldname': edge, 'weight': H.get_hyperedge_weight(edge)})
            eid = H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge),weight =H.get_hyperedge_weight(edge))
            #print('tail of edge {}'.format(eid))
            #if len(H.get_hyperedge_tail(edge)) < 20:
                #for v in H.get_hyperedge_tail(edge):
                    #print(v[-4:])
            #print('head of edge {}'.format(eid))
            #if len(H.get_hyperedge_head(edge)) < 20:
                #for v in H.get_hyperedge_head(edge):
                    #print(v[-4:])
        #print('done adding new edges')
        #print('target',target, H2.get_backward_star(target))
        H2.add_node('SUPERSOURCE',{'label': 'SUPERSOURCE'}) 
        for edge in H2.get_forward_star('SUPERSOURCE'):
            #print(edge,H2.get_
            forwardstarlist = []
            for v in H2.get_hyperedge_head(edge):
                #print(v,H2.get_forward_star(v))
                if len(H2.get_forward_star(v)) > 0:
                    forwardstarlist.append(v)
            #print('old list',H2.get_hyperedge_head(edge))
            H2.remove_hyperedge(edge)
            #print('new list',forwardstarlist)
            if len(forwardstarlist) > 0:
                H2.add_hyperedge(['SUPERSOURCE'],forwardstarlist,weight=0)
            else:
                H2.add_hyperedge(['SUPERSOURCE'],[],weight=0)

        H = H2
    else:
        print('already have smallest edgeset')
    #print('done removing nodes with no outedges from the source')
    #NOTE: This was removed for the enumeration code and needs to be added back in 
    #for v in H.get_node_set():
        ##need to remap the edges because just calling remove_node(v) removes also all the hyperedges v is in
        #if len(H.get_forward_star(v)) == 0 and v != 'SUPERTARGET':
            ##find edges that need to be replaced
            ##H.trim_node(v)
            #backedges = H.get_backward_star(v)
            #for e in backedges:
                #tail = H.get_hyperedge_tail(e)
                #head = H.get_hyperedge_head(e)
                #head.remove(v)
                #w = H.get_hyperedge_weight(e)
                ##oldattrs = H.get_hyperedge_attributes(e)
                #H.remove_hyperedge(e)
                ##H.add_hyperedge(tail,head,oldattrs)
                #H.add_hyperedge(tail,head,weight=w)
            #H.remove_node(v)
    ##for v in H.get_node_set():
        ##if len(H.get_forward_star(v)) == 0 and v != 'SUPERTARGET':
            ##print('There has been a huge problem with vertex',v)
    ##print('done removing nodes with no outedges')
    #NOTE: End of taking out
    #printremappededges(H)
    #print('building up reachable edges')
    reachableedges = []
    #H2.add_node('SUPERSOURCE',{'label': 'SUPERSOURCE'}) 
    for edge in H.hyperedge_id_iterator():
        reachableedges.append(edge)

    #initialize edgedict
    #for e in reachableedges:
    for e in H.hyperedge_id_iterator():
        edgedict[e] = {'isremoved': False, 'bestinedges': [], 'candidateinedges': [], 'tailcount': len(H.get_hyperedge_tail(e))}
    #print('initialized edge dict')
    #print(edgedict)

    #initialize taildistancelist
    taildistancelist = {}
    for e in reachableedges:
        taildistancelist[e] = 'inf'

    #initialize reachableedgecounter
    reachableedgecounter = len(reachableedges)
    print(('reachableedgecounter: {0}'.format(reachableedgecounter)))
    #print(('{0} nodes in hypergraph'.format(len(H.get_node_set()))))

    #print('\n\n\n\n\n\n\n',[H.get_hyperedge_tail(e) for e in H.get_backward_star(6145)])
    #initialize heap
    heap = []
    entry_finder = {}               # mapping of tasks to entries, this and the line below are strictly for the heap
    counter = itertools.count()     # unique sequence count
    for e in H.get_forward_star(source):
        if e in reachableedges:
            add_node(heap,e,weight(H,[e]),counter,entry_finder)
    #print('heap initialized')
    #print(heap)

    #initialize reachedtable
    reachedtable = {}
    for e in H.hyperedge_id_iterator():
        reachedtable[e] = []

    return reachableedges,edgedict,taildistancelist,heap,reachableedgecounter,reachedtable, entry_finder, counter, H

def printremappededges(H):
    print('starting to remap edges')
    nodedict = {}
    nodedict['SUPERSOURCE'] = 'SUPERSOURCE'
    nodedict['SUPERTARGET'] = 'SUPERTARGET'
    #ofile = open('out','r').readlines()
    ofile = json.load(open('physical_entities.json','r'))
    #for i in range(0,len(ofile),2):
        #nodedict[ofile[i].strip()] = ofile[i+1].strip()
    for entity in ofile:
        H.add_node(entity['uri'],{'label': entity['name']})
        nodedict[entity['uri']] = entity['name']
    for v in H.node_iterator():
        try:
            H.get_node_attribute(v,'label')
        except:
            H.add_node(v,{'label': v})
    for e in H.hyperedge_id_iterator():
        #print(e)
        tail = []

        for v in H.get_hyperedge_tail(e):
            if v in nodedict:
                tail.append(nodedict[v])
            else:
                nodedict[v] = v
                tail.append(v)
            tail.append(H.get_node_attribute(v,'label'))
        #print('tail',tail)
        head = []
        for v in H.get_hyperedge_head(e):
            if v in nodedict:
                head.append(nodedict[v])
            else:
                nodedict[v] = v
                head.append(v)
            head.append(H.get_node_attribute(v,'label'))
        #print('head',head)
    return nodedict

def findreachableandbackrecoverable(H,source,target):
    '''
    Finds and returns the set of reachable and backwards-recoverable edges
    '''
    reachable = findreachable(H,source)
    print(H.get_backward_star(target))
    targete = [e for e in H.get_backward_star(target) if e in reachable]
    print('targete',targete)
    #print('\n\n\n\n\n\n\n',[H.get_hyperedge_tail(e) for e in H.get_backward_star(6145)])
    #print('reachable set')
    #print(reachable)
    print(('{0} reachable edges'.format(len(reachable))))

    recoverable = findrecoverable(H,target)
    #print('recoverable set')
    #print(recoverable)
    print(('{0} recoverable edges'.format(len(recoverable))))
    
    reachrecover = list(set(reachable) & set(recoverable))
    return reachrecover

def findreachable(H,source):
    '''
    Finds the set of reachable edges from the source
    '''
    _,nodes,edges,___ = b_visit(H,source)
    if 'SUPERTARGET' in nodes:
        print('supertarget here',nodes['SUPERTARGET'])
        for edge in H.get_backward_star('SUPERTARGET'):
            for node in H.get_hyperedge_tail(edge):
                if nodes[node] == None:
                    print(node)
    else:
        print('supertarget not here')
    reachable = []
    for e in edges:
        if edges[e] != None:
            reachable.append(e)
    
    return reachable

def findrecoverable(H,target):
    '''
    Finds the set of backwards-recoverable edges from the sink
    '''
    sinkedges = []
    for e in H.get_backward_star(target):
        sinkedges.append(e)

    recoverable = set(sinkedges)
    print(('there are {0} sink edges'.format(len(sinkedges))))
    for e in sinkedges:
        print((e,H.get_hyperedge_tail(e),'head',H.get_hyperedge_head(e),[H.get_backward_star(v) for v in H.get_hyperedge_tail(e)]))
        _,R,__ = recover(H,e,'all',{})
        #print(('from {0} I recovered {1} edges'.format(e,len(R))))
        for f in R:
            recoverable.add(f)
    return recoverable
        

def reachability_from_edges(H,ones,source):
    #takes the solution from the previous iteration of the ILP and sees what nodes are reachable from those edges.  Make a subhypergraph, and then call b_visit on it.
    #returns a list of nodes, S
    print('_____________rfe_______')
    print(ones)
    nodeset = set()
    H2 = DirectedHypergraph()
    H2.add_nodes(H.get_node_set())
    for edge in ones:
        H2.add_hyperedge(H.get_hyperedge_tail(una(edge)),H.get_hyperedge_head(una(edge)))
    
    nodeset = b_visit(H2,source)
    print('nodeset size:',len(nodeset[0]), nodeset[0])
    return nodeset


def add_node(pq,node, priority, counter, entry_finder):
    'Add a new task or update the priority of an existing task'
    if node in entry_finder:
        remove_node(pq,node,entry_finder)
    count = next(counter)
    entry = [priority, count, node]
    entry_finder[node] = entry
    heapq.heappush(pq, entry)

def remove_node(pq,node,entry_finder):
    'Mark an existing task as REMOVED.  Raise KeyError if not found.'
    REMOVED = '<removed-node>'      # placeholder for a removed node
    entry = entry_finder.pop(node)
    entry[-1] = REMOVED

def pop_node(pq,entry_finder):
    'Remove and return the lowest priority task. Raise KeyError if empty.'
    REMOVED = '<removed-node>'      # placeholder for a removed node
    while pq:
        priority, count, node = heapq.heappop(pq)
        if node != REMOVED:
            del entry_finder[node]
            return node,priority
    raise KeyError('pop from an empty priority queue')

def emptypq(pq,entry_finder):
    REMOVED = '<removed-node>'      # placeholder for a removed node
    while pq:
        priority, count, node = heapq.heappop(pq)
        if node != REMOVED:
            heapq.heappush(pq,[priority, count, node])
            return False
    return True
    
        
def sort_edges_in_path(H,path,source):
    '''
    Here we want to put the hyperedges in a path into an order that works for the initial
    set of inequalities.  To do this, we start with just s in S, then add to S any head of 
    a hyperedge that is reachable with just the edges in S, and repeat
    '''
    orderedpath = []
    S = set()
    S.add(source)
    pathlen = len(path)
    while path:
        for edge in path:
            for node in H.get_hyperedge_tail(edge):
                if node not in S:
                    break
            else:
                orderedpath.append(edge)
                path.remove(edge)
                for headnode in H.get_hyperedge_head(edge):
                    S.add(headnode)
    return orderedpath
    
def visit_set(H,x,ilp,source):
    y = ilp.solution.get_values(x)
    hedges = [x[i].replace('a_','') for i in range(len(x)) if y[i]==1]
    HSUB = DirectedHypergraph()
    for hedge in hedges:
        HSUB.add_hyperedge(H.get_hyperedge_tail(hedge),H.get_hyperedge_head(hedge),weight=1)
    nodes,ignore,ignore = visit(HSUB,source)
    return nodes,hedges

