from xml.dom import minidom
import sys
import cplex
import itertools
import heapq
import time
import pickle as pkl
from multiprocessing import Process, Value, Array

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit,visit
from halp.utilities import directed_graph_transformations

from pathheuristic import trim,weight
from cutfinder import *
#need to get the accumulation and target production constraints, plus the edge constraints that tie real value to indicator variables

EPSILON=0.0001
CONSTANT=1000000

def make_shortest_cyclic_hyperpath_ilp(H,source,target,outfile):
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')
    c = 0
    #out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
    #out,c = writeConstraint_IfHnodeThenBackwardStar(H,source,out,c)
    #out,c = writeConstraint_FixValues(H,{target:1},out,c)
    #out,c = writeConstraint_OrderVariables(H,out,c)

    #THIS IS WHAT MAKES IT THE ILP

    #out = writeBinaryBounds(H,out)

    out.write('End\n')
    out.close()


def make_shortest_acyclic_hyperpath_ilp(H,source,target,outfile):
    for node in H.node_iterator():
        if type(node) is not int:
            print(node)
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')
    c = 0
    out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
    out,c = writeConstraint_IfHnodeThenBackwardStar(H,source,out,c)
    out,c = writeConstraint_FixValues(H,{target:1},out,c)
    out,c = writeConstraint_OrderVariables(H,out,c)
    out = writeBinaryBounds(H,out)
    out.write('End\n')
    out.close()

    print(('Wrote to %s' % (outfile)))
    return

def make_shortest_test(H,source,target,outfile,high_penalty_sources):
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')
    c = 0
    out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
    out,c = writeConstraint_IfHnodeThenBackwardStar(H,source,out,c)
    out,c = writeConstraint_FixValues(H,{target:1},out,c)
    out,c,G = writeConstraint_Flow(H,source,target,out,c,high_penalty_sources=high_penalty_sources)
    out = writeBinaryBounds(H,out,write_f=G)
    out.write('End\n')
    out.close()

    print(('Wrote to %s' % (outfile)))
    return

def make_shortest_hyperpath_ilp_predecessors(H,source,target,outfile):
    '''
    BROKEN - doesn't work, but keep the code anyway.
    '''
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')
    c = 0
    out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
    out,c = writeConstraint_IfHedgeThenPredecessor(H,out,c)
    out,c = writeConstraint_IfPredecessorThenHedge(H,out,c)
    out,c = writeConstraint_IfHnodeThenPredecessor(H,[source],out,c)
    out,c = writeConstraint_fixTails([target],out,c)
    out = writeBounds(H,out,targets=[target])

    out.write('End\n')
    out.close()

    print(('Wrote to %s' % (outfile)))
    return

def make_shortest_hyperpath_ilp_simple(H,source,target,outfile):
    out = open(outfile,'w')

    ## write objective
    out = writeObjective(H,out)

    ## write constraints
    out.write('Subject To\n')
    c = 0
    out,c = writeConstraint_IfHedgeThenIncidents(H,out,c)
    out,c = writeConstraint_IfHnodeThenBackwardStar(H,source,out,c)
    out,c = writeConstraints_SimpleWalks(H,out,c,source,target)
    out,c = writeConstraint_FixValues(H,{target:1},out,c)
    out = writeBinaryBounds(H,out,write_q=True)

    out.write('End\n')
    out.close()

    print(('Wrote to %s' % (outfile)))
    return

def writeObjective(H,out,minimize=True):
    if minimize:
        out.write('Minimize\n')
    else:
        out.write('Maximize\n')

    for hedge in H.hyperedge_id_iterator():
        out.write(' + %d %s' % (H.get_hyperedge_attribute(hedge,'weight'),a(hedge)))
    out.write('\n')
    return out

def writeConstraint_IfHedgeThenIncidents(H,out,c):
    '''
    Writes the constraint that if hyperedge e is in the 
    solution, then all hypernodes incident to e must also
    be in the solution.  That is, for $I(e) = H(e) \cup T(e)$,
    \sum_{u I(e)} a_u >= |I(e)| a_e 
    \sum_{u I(e)} a_u - |I(e)| a_e >= 0

    '''
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        incident_set = set(H.get_hyperedge_head(hedge)).union(H.get_hyperedge_tail(hedge))
        out.write('c%d_if_hedge_then_incidents: ' % (c))
        for hnode in incident_set: # for all u \in T(e) \cup H(e)
            #out.write('+ %s ' % a(hnode))
            out.write('+ %s ' % a(hnode))
        out.write(' - %d %s >= 0\n' % (len(incident_set),a(hedge)))
        c+=1

    return out,c

def writeConstraint_IfHnodeThenBackwardStar(H,s,out,c):
    '''
    Writes the constraint that if hypernode u is in the solution,
    then there must be at least one hyperedge in the backwards star
    that is also in the solution. This holds for all hypernodes
    except for s.
    \sum_{e \in BS(v)} a_e >= a_v
    \sum_{e \in BS(v)} a_e - a_v >= 0
    '''
    for hnode in H.node_iterator():
        if hnode == s:
            continue
        out.write('c%d_if_hnode_then_backwardstar: ' % (c))
        for hedge in H.get_backward_star(hnode):
            out.write('+ %s ' % (a(hedge)))
        out.write('- %s >= 0\n' % (a(hnode)))
        c+=1

    return out,c

def writeConstraint_FixValues(H,tofix,out,c):
    '''
    Writes the constraint that fixes variables in the tofix dictionary.
    The variable tofix is a dictionary of {hypernode: <0 or 1>}.
    '''
    for t in tofix:
        out.write('c%d_fixed: %s = %d\n' % (c,a(t),tofix[t]))  
        c+=1 #increment constraint counter
    return out,c

'''
def writeConstraint_AnyTarget(targets,out,c):
    out.write('c%d_anytarget: ' % (c))
    for t in targets:
        out.write('+ %s' % (a(t)))
    out.write(' >= 1\n')  
    c+=1 #increment constraint counter
    return out,c
'''
def writeConstraint_HyperedgeOrderVariables(H,out,c,sources):
    ''' 
    Writes the constraint that the at least one incoming order
    variable must be smaller than at least one outgoing order variable.
    '''

    '''
    Order Bounds: a_v >= o_v >= 0
    o_v >= 0 and a_v - o_v >= 0
    '''
    for hnode in H.node_iterator():
        out.write('c%d_order_bounds: %s >= 0\n' % (c,o(hnode)))
        c+=1 

        out.write('c%d_order_bounds: %s - %s >= 0\n' % (c,a(hnode),o(hnode)))
        c+=1 #increment constraint counter

    '''
    Order Bounds: a_e >= o_e >= 0
    e_e >= 0 and a_e - o_e >= 0
    '''
    for hedge in H.hyperedge_id_iterator():
        out.write('c%d_order_bounds: %s >= 0\n' % (c,o(hedge)))
        c+=1 

        out.write('c%d_order_bounds: %s - %s >= 0\n' % (c,a(hedge),o(hedge)))
        c+=1

    '''
    Order value of hypernode is MAX of order value on all OUTGOING hyperedges.
    o_v >= o_e forall v \in V, e \in FS(v)
    o_v - o_e >= 0 forall v \in V, e \in FS(v)
    '''
    for v in H.node_iterator():
        for e in H.get_forward_star(v):
                out.write('c%d_hyperenode_order: %s - %s >= 0\n' % (c,o(v),o(e)))
                c+=1 
    '''
    Finally, orrder value for INCOMING hyperedges is smaller than order value 
    of any node in the HEAD.
    o_e < o_v + C (1-a_e) forall e \in E, v \in H(e)
    o_e <= o_v + C (1-a_e) - epsilon forall e \in E, v \in H(e)
    o_e - o_v + C a_e <= C - epsilon forall e \in E, v \in H(e)
    '''
    for e in H.hyperedge_id_iterator():
        for v in H.get_hyperedge_head(e):
            out.write('c%d_hyperedge_order: %s - %s + %d %s <= %f\n' % (c,o(e),o(v),CONSTANT,a(e),CONSTANT-EPSILON))
            c+=1

    return out,c

def writeConstraint_OrderVariables(H,out,c):
    ''' 
    Writes the constraint that the order variables
    must be smaller in the tail than the order
    variables in the head.

    Order Bounds: a_v >= o_v >= 0
    o_v >= 0 and a_v - o_v >= 0
    '''
    for hnode in H.node_iterator():
        out.write('c%d_order_bounds: %s >= 0\n' % (c,o(hnode)))
        c+=1 

        out.write('c%d_order_bounds: %s - %s >= 0\n' % (c,a(hnode),o(hnode)))
        c+=1 #increment constraint counter

    '''
    Order Constraints:
    o_u <= o_v - \epsilon + C (1-a_e) forall u,v in (T(e),H(e)) forall e \in E 
    o_u <= o_v - \epsilon + C (1-a_e) 
    o_u - o_v + C a_e <= C - \epsilon
    '''
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        for u in H.get_hyperedge_tail(hedge):
            for v in H.get_hyperedge_head(hedge):
                out.write('c%d_order_constraint: %s - %s + %d %s <= %f\n' % \
                              (c,o(u),o(v),CONSTANT,a(hedge),CONSTANT-EPSILON))
                c+=1

    return out,c

def writeConstraint_Flow(H,source,target,out,c,high_penalty_sources=set()):
    '''
    Converts H to a graph G and imposes a flow constraint 
    '''
    G = directed_graph_transformations.to_graph_decomposition(H)
    graph_edges = [(G.get_hyperedge_tail(hedge)[0],G.get_hyperedge_head(hedge)[0]) for hedge in G.get_hyperedge_id_set()]
    print(('%d graph edges' % (len(graph_edges))))
    '''
    Require that exactly one outgoing edge from s has an f value of 1.
    \sum_{(u,v): u = s} f_{uv} = 1
    '''
    print ('1...')
    out.write('c%d_outgoing_flow: ' % (c))
    for e in [e2 for e2 in graph_edges if e2[0]==source and e2[1] not in high_penalty_sources]:
        out.write(' + %s' % (f(e[0],e[1])))
    out.write(' = 1\n')
    c+=1

    '''
    Require that exactly one incoming edge to t has an f value of 1.
    \sum_{(u,v): v = t} f_{uv} = 1
    '''
    print ('2...')
    out.write('c%d_incoming_flow: '  %(c))
    for e in [e2 for e2 in graph_edges if e2[1]==target]:
        out.write(' + %s' % (f(e[0],e[1])))
    out.write(' = 1\n')
    c+=1

    '''
    Require that each non-s/t node is balanced.
    \sum_{e \in BS(v)} f_e  = \sum_{e' \in FS(v)} f_e' forall v \in V \setminus {s,t}
    \sum_{e \in BS(v)} f_e - \sum_{e' \in FS(v)} f_e' = 0 forall v \in V \setminus {s,t}
    '''
    print ('2.5...')
    #BS = {u:set() for u in G.node_iterator()}
    #FS = {u:set() for u in G.node_iterator()}
    BS = {}
    FS = {}
    for u in G.node_iterator():
        BS[u] = set()
        FS[u] = set()
    for u,v in graph_edges:
        FS[u].add(v)
        BS[v].add(u)
    print ('3...')
    for v in G.node_iterator():
        if v == target or v == source:
            continue # skip constraint 
        out.write('c%d_flow_balance_%s: ' % (c,v))
        for u in BS[v]:
            out.write(' + %s' %(f(u,v)))
            for w in FS[v]:
                out.write(' - %s' %(f(v,w)))
        out.write(' = 0\n')
        c+=1

    '''
    Require that if f is 1, then alpha is 1.
    f_{uv} <= \sum_{e: u \in T(e), v \in H(e)} a_e forall u,v in Graph edges
    \sum_{e: u \in T(e), v \in H(e)} a_e - f_{uv} >= 0 forall u,v in Graph edges
    '''
    print ('3.5...')
    hedge_mapper = {}
    for e in graph_edges:
        hedge_mapper[e] = set()
    #hedge_mapper = {e:set() for e in graph_edges}
    for hedge in H.hyperedge_id_iterator():
        for u in H.get_hyperedge_tail(hedge):
            for v in H.get_hyperedge_head(hedge):
                hedge_mapper[(u,v)].add(hedge)
    print ('4...')
    for u,v in graph_edges:
        out.write('c%d_flow_hedge_match:' % (c))
        for hedge in hedge_mapper[(u,v)]:
            out.write(' + %s' % (a(hedge)))
        out.write(' - %s >= 0\n' % (f(u,v)))
        c+=1
    return out,c,G

def writeBinaryBounds(H,out,write_q=False,write_f=None):
    '''
    Specify all the alpha variables as binary.
    '''
    out.write('Binary\n')
    for hnode in H.node_iterator():
        out.write(' %s\n' % (a(hnode)))
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        out.write(' %s\n' % (a(hedge)))
        if write_q:
            out.write(' %s\n' % (q(hedge)))
    if write_f:
        ## write_f is a GRAPH
        edges = [(write_f.get_hyperedge_tail(e)[0],write_f.get_hyperedge_head(e)[0]) for e in write_f.get_hyperedge_id_set()]
        for e in edges:
            out.write(' %s\n' % (f(e[0],e[1])))
    return out

def writeBinaryBoundsCycles(H,out):
    '''
    Specify all the edge variables as binary.
    '''
    out.write('Binary\n')
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        out.write(' %s\n' % (a(hedge)))
    return out

####################
def writeConstraint_IfHedgeThenPredecessor(H,out,c):
    '''
    Writes the constraint that if hyperege e is in the solution,
    then it must act as at least one predecessor for a hypernode
    in the head of e. That is,
    a_e \leq \sum_{u \in H(e)} \pi_{e,u}         for all e \in E
    a_e - \sum_{u \in H(e)} \pi_{e,u} \leq 0     for all e \in E
    '''
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        out.write('c%d_if_hedge_then_predecessor: %s' % (c,a(hedge)))
        for hnode in H.get_hyperedge_head(hedge):
            out.write(' - %s' % (pi(hnode,hedge)))
        out.write(' <= 0\n')
        c+=1

    return out,c

def writeConstraint_IfPredecessorThenHedge(H,out,c):
    '''
    Writes the constraint that if a predecessor variable is in
    the solution, then the hyperege e is in the solution.  That is,
    \pi_{e,u} \leq a_e             for all e \in E, u \in H(e)
    \pi_{e,u} - a_e \leq 0        for all e \in E, u \in H(e)
    '''
    for hedge in H.hyperedge_id_iterator(): # for all e \in E
        for hnode in H.get_hyperedge_head(hedge): #\forall u \in H(e)
            out.write('c%d_if_predecessor_then_hedge: %s - %s <= 0\n' % (c,pi(hnode,hedge),a(hedge)))
            c+=1

    return out,c

def writeConstraint_IfHnodeThenPredecessor(H,sources,out,c):
    '''
    Writes the constraint that if hypernode u is in the solution,
    then there must be exactly one predecessor variable specified.
    This only holds true for hypernodes that are not in S.
    a_u = \sum_{e:u \in H(e)} \pi_{e,u}         for all u \in U \setminus S
    a_u - \sum_{e:u \in H(e)} \pi_{e,u} = 0        for all u \in U \setminus S
    '''
    for hnode in H.node_iterator():
        if hnode in sources:
            continue
        out.write('c%d_if_hnode_then_predecessor: %s' % (c,a(hnode)))
        for hedge in H.get_backward_star(hnode):
            out.write(' - %s' % (pi(hnode,hedge)))
        out.write(' = 0\n')
        c+=1

    return out,c

def writeConstraint_fixTails(targets,out,c):
    '''
    Sets alpha values of targets to 1.
    a_u = 1         for all u \in T
    '''
    for hnode in targets:
        out.write('c%d_fix_tails: %s = 1\n' % (c,a(hnode)))
        c+=1
    return out,c

def writeBounds(H,out,targets=None):
    out.write('Binary\n')
    for hnode in H.node_iterator():
        out.write(a(hnode)+'\n')
    for hedge in H.hyperedge_id_iterator():
        out.write(a(hedge)+'\n')
        if targets:
            for t in targets:
                out.write(pi(t,hedge)+'\n')
        else:
            for hnode in H.get_hyperedge_head(hedge):
                out.write(pi(hnode,hedge)+'\n')
    return out

####################
def writeConstraints_SimpleWalks(H,out,c,source,target):
    '''
    Writes three types of constraints with q variables to get simple walks.
    '''

    '''
    Require that exactly one outgoing hyperedge from has a q value of 1.
    \sum_{e:s \in T(e)} q_e = 1
    '''
    out.write('c%d_outgoing_walk: ' % (c))
    for hedge in H.get_forward_star(source):
        out.write(' + %s' % (q(hedge)))
    out.write(' = 1\n')
    c+=1

    '''
    Require that exactly one incoming hyperedge to t has a q value of 1.
    \sum_{e:t \in H(e)} q_e = 1
    TODO: deal with multiple targets - does this still work?
    '''
    out.write('c%d_incoming_walk: '  %(c))
    for hedge in H.get_backward_star(target):
        out.write(' + %s'  % (q(hedge)))
    out.write(' = 1\n')
    c+=1

    '''
    Require that the heads of each hyperedge is balanced: here "balanced"
    means that q is <= the summ of the q values for all outgoing hyperedges
    from the heads of e.
    q_e <= \sum_{v:v \in H(e)} \sum_{e':v \in T(e'); e' != e} q_e'     forall e \in E \setminus \{e:t \in H(e)}
    \sum_{v:v \in H(e)} \sum_{e':v \in T(e'); e' != e} q_e' - q_e >= 0     forall e \in E \setminus \{e:t \in H(e)\}
    '''
    for hedge in H.hyperedge_id_iterator():
        if target in H.get_hyperedge_head(hedge):
            continue # skip constraint for hyperedges coming into target.
        out.write('c%d_simple_balance_outgoing_%s: ' % (c,hedge))
        for hnode in H.get_hyperedge_head(hedge):
            for next_hedge in H.get_forward_star(hnode):
                if next_hedge != hedge:
                    out.write(' + %s' %(q(next_hedge)))
        out.write(' - %s >= 0\n' % (q(hedge)))
        c+=1

    '''
    Require that the tails of each hyperedge is balanced: here "balanced"
    means that q is <= the sum of the q values for all incoming hyperedges
    to the tails of e.
    q_e <= \sum_{u:u \in T(e)} \sum_{e':u \in H(e'); e' != e} q_e'     forall e \in E \setminus \{e:s \in T(e)}
    \sum_{u:u \in T(e)} \sum_{e':u \in H(e'); e' != e} q_e' - q_e >= 0     forall e \in E \setminus \{e:s \in T(e)\}
    '''
    for hedge in H.hyperedge_id_iterator():
        if source in H.get_hyperedge_tail(hedge):
            continue # skip constraint for hyperedges coming out of source.
        out.write('c%d_simple_balance_incoming_%s: ' % (c,hedge))
        for hnode in H.get_hyperedge_tail(hedge):
            for prev_hedge in H.get_backward_star(hnode):
                if prev_hedge != hedge:
                    out.write(' + %s' %(q(prev_hedge)))
        out.write(' - %s >= 0\n' % (q(hedge)))
        c+=1

    '''
    Require that if q is 1, then alpha is 1.
    q_e <= a_e             forall e \in E
    q_e - a_e <= 0      forall e \in E
    '''
    for hedge in H.hyperedge_id_iterator():
        out.write('c%d_if_q_then_a: %s - %s <= 0\n' % (c,q(hedge),a(hedge)))
        c+=1

    return out,c

def t_in_S(S,target):
    #here is some code about checking if the sink is in S
    #returns a boolean
    for node in S:
        if node == target:
            return True
    return False

def get_new_constraint(S,H,headcut=None):
    #this takes the S,T cut and the set of edges and finds all edges that cross the cut
    #this one is a little harder to think through, but this is just the exponential version.  It should be fine for our hypergraphs, but does not work for any hypergraph.
    #returns a list of edges that cross the s,t-cut
    #print( '_____________S__________',S)
    crossedges = []
    addedge = False
    nodesintail = True
    for hedge in H.hyperedge_id_iterator():
        addedge = False
        nodesintail = True
        for tailnode in H.get_hyperedge_tail(hedge):
            if tailnode not in S:
                nodesintail = False
        #if it got here, that means that all the nodes in the tail are in S
        if nodesintail == True:
            for headnode in H.get_hyperedge_head(hedge):
                if headnode not in S:
                    #this means that this edge crosses the cut
                    addedge = True
            if addedge == True:
                crossedges.append(hedge)
    return crossedges


def sourcecut_from_headcut(H,activeedges,headcut):
    nodeset = set()
    H2 = DirectedHypergraph()
    H2.add_nodes(H.get_node_set())
    for edge in activeedges:
        H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge))

    source_side_vertices = set()
    for e in activeedges:
        for v in H.get_hyperedge_tail(e):
            if v in headcut:
                source_side_vertices.add(v)

    newsource = 'NEWTEMPSOURCE'
    H2.add_hyperedge(set([newsource]),source_side_vertices)
            

    nodeset,t1,t2,t3 = b_visit(H2,newsource)
    newheadcut = set(headcut)
    for v in nodeset:
        newheadcut.add(v)
    return tuple(newheadcut)

def find_new_active_crossing_edges(H,activeedges,crossing_active_edges,headcut,v):

    new_active_crossing_count = 0
    new_active_crossing_edges = []
    #We now want to count what new edges will be crossing if we were to move v to the sink side
    for f in H.get_backward_star(v):
        if f in activeedges and f not in crossing_active_edges:
            #determine if v moving to the sink would make f newly crossing
            tail_in_source = True
            for w in H.get_hyperedge_tail(f):
                if w not in headcut or w == v:
                    tail_in_source = False
                    break
            if tail_in_source == True:
                #If the tail is all on the source side, moving v (which is in e's head) would make it now crossing
                if f not in new_active_crossing_edges:
                    new_active_crossing_count += 1
                    new_active_crossing_edges.append(f)
    return new_active_crossing_count,new_active_crossing_edges
                            

def sinkcut_from_headcut(H,activeedges,headcut):
    #first we need to determine the initial crossing hyperedges
    crossing_active_edges = set()

    for e in activeedges:
        tail_in_source = True
        for v in H.get_hyperedge_tail(e):
            if v not in headcut:
                tail_in_source = False
                break
        if tail_in_source == True:
            for v in H.get_hyperedge_head(e):
                if v not in headcut:
                    crossing_active_edges.add(e)
                    break

    newheadcut = set(headcut)
    while len(crossing_active_edges) > 0:
        #print('another while loop iteration')
        min_hyperedge_cross_from_vertex = None
        min_hyperedge_cross_vertex = None
        min_new_active_crossing_edges = []
        for e in crossing_active_edges:
            #Check each tail vertex of the crossing_active_edges to see if moving them over gives the fewest new crossing edges
            for v in H.get_hyperedge_tail(e):
                if v in newheadcut:
                    new_active_crossing_count, new_active_crossing_edges = find_new_active_crossing_edges(H,activeedges,crossing_active_edges,newheadcut,v)
                    #Now we want to see if v moves the fewest active edges so far
                    if min_hyperedge_cross_from_vertex == None or new_active_crossing_count < min_hyperedge_cross_from_vertex:
                        min_hyperedge_cross_from_vertex = new_active_crossing_count
                        min_hyperedge_cross_vertex = v
                        min_new_active_crossing_edges = new_active_crossing_edges

        if min_hyperedge_cross_vertex != None:
            newheadcut.remove(min_hyperedge_cross_vertex)

        #find the hyperedges that no longer cross the cut and then remove them
        former_crossing_edges = []
        for e in crossing_active_edges:
            if min_hyperedge_cross_vertex in H.get_hyperedge_tail(e):
                former_crossing_edges.append(e)

        #print('len former crossing cuts is {}'.format(len(former_crossing_edges)))
        for e in former_crossing_edges:
            crossing_active_edges.remove(e)

        for e in min_new_active_crossing_edges:
            crossing_active_edges.add(e)

    return tuple(newheadcut)


def reachability_from_edges(H,ones,source):
    #takes the solution from the previous iteration of the ILP and sees what nodes are reachable from those edges.  Make a subhypergraph, and then call b_visit on it.
    #returns a list of nodes, S
    #print( '_____________rfe_______')
    #print( ones)
    nodeset = set()
    H2 = DirectedHypergraph()
    H2.add_nodes(H.get_node_set())
    for edge in ones:
        H2.add_hyperedge(H.get_hyperedge_tail(una(edge)),H.get_hyperedge_head(una(edge)))
    
    nodeset = b_visit(H2,source)
    #print ('nodeset size:',len(nodeset[0]), nodeset[0])
    return nodeset

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
        print(( curnode,curpriority))
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
    
def inequalities_from_path(H,path,source,target):
    I = []
    C = set()
    C.add(source)
    visited_nodes = {}
    path = sort_edges_in_path(H,path,source)
    usededges = []
    while path and target not in C:
        print(( 'major cut', C))
        inequalraw,inequal,nodes_to_add,minoredges,currusededges = major_cut(H,path,C,visited_nodes)
        I.append(inequal)
        for minoredge in minoredges:
            print( minoredge)
            minorinequal,Cprime = minor_cut(H,minoredge,C)
            print(( 'minor cut', Cprime))
            I.append(minorinequal)

        C.update(nodes_to_add)
        removefrompath = []
        newusededges = []
        usededges = usededges + currusededges
        for i in range(len(usededges)):
            edge,count = usededges[i]
            if count == 1:
               removefrompath.append(edge) 
            else:
                newusededges.append((edge,count - 1))
        usededges = newusededges
                
        path = [x for x in path if x not in removefrompath]
    return I

def major_cut(H,path,C,visited_nodes):
    '''
    Takes an edge and gives the major cut as an inequality and the nodes that will be added
    to C after the major cut is added to I.  Also outputs the set of edges crossing the major
    cut that aren't the edge in the path to make it easy to make minor edges
    '''
    #find the edges that cross the cut
    print( C)
    verticestoadd = []
    crossedges = []
    minoredges = []
    usedmajors = []
    possibleminors = []
    for edge in H.hyperedge_id_iterator():
        for node in H.get_hyperedge_tail(edge):
            if node not in C:
                break
        else:
            for node in H.get_hyperedge_head(edge):
                if node not in C:

                    crossedges.append(edge)
                    break
    count = 0
    for edge in crossedges:
        if edge in path:
            count += 1
            possibleminors.append(edge)
        else:
            minoredges.append(edge)
    for edge in crossedges:
        for node in H.get_hyperedge_head(edge):
            if edge in path:
                usedmajors.append((edge,count))
                if node in list(visited_nodes.keys()):
                    visited_nodes[node] -= 1
                else:
                    visited_nodes[node] = count
            else:
                visited_nodes[node] = 1
    for node in list(visited_nodes.keys()):
        if visited_nodes[node] == 1:
            verticestoadd.append(node)
    inequal = make_constraint_from_edges(crossedges)
    if count > 1:
        minoredges = minoredges + possibleminors
    return crossedges,inequal,verticestoadd,minoredges,usedmajors

def get_addable_cuts(cuts,realvaluedcuts=False):
    '''

    '''
    constraints = []
    for cut in cuts:
        #print(cut)
        constraint = make_constraint_from_edges(cut,realvaluedcuts=realvaluedcuts)
        #print(constraint)
        constraints.append(constraint)
    return constraints

def make_constraint_from_edges(crossedges,realvaluedcuts=False):
    '''

    '''

    #addconstraint to the ILP
    crossedges = set(crossedges)
    if realvaluedcuts == True:
        ones = [fedge(c) for c in crossedges]
    else:
        ones = [a(c) for c in crossedges]
    val = [1.0]*len(ones)
    eq = cplex.SparsePair(ind=ones,val=val)
    return eq

            
            
            
    #modify the queue to see if any nodes have expired (for now I'm just implementing this as a list of tuples) 


def minor_cut(H,minoredge,C):
    '''
    Gives a minor cut inequality taking in the major cut inequality and the minor edge
    Just searches for the edges that need to be added when the minor cut is made to not cross the cut.
    '''
    Cprime = C.union(set(H.get_hyperedge_head(minoredge)))

    crossedges = []
    for edge in H.hyperedge_id_iterator():
        for node in H.get_hyperedge_tail(edge):
            if node not in Cprime:
                break
        else:
            for node in H.get_hyperedge_head(edge):
                if node not in Cprime:
                    crossedges.append(edge)
    inequal = make_constraint_from_edges(crossedges)
    return inequal,Cprime
        
        
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
    
def write_and_print_vars(ilp,H,nodeset,numsolsfound,verbose,outprefix,ilpprint=False):
    '''

    '''
    print(('status',ilp.solution.get_status()))
    print(('feasible',ilp.solution.is_primal_feasible()))
    print(('objective value',ilp.solution.get_objective_value()))
    print(('-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
    
    print ('Solution Found.')
    #extract the edges from the solution
    objective = ilp.solution.get_objective_value()
    ilp.solution.write('%s-%d.sol' % (outprefix,numsolsfound))
    variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
    ones = [var for var in variables if variables[var] > 0 and var[0:3]=='a_e' and var.split('_')[1] not in nodeset]
    if ilpprint == True:
        unaones = [una(e) for e in ones]
        for e in unaones:
            print((e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e)))
    print((ones,'ones'))

def ILP_from_LP(lp,H,nodeset):
    '''

    '''
    #print('-'*30,'ILP')
    for e in H.hyperedge_id_iterator():
        lp.variables.set_types(a(e),lp.variables.type.binary)
    lp.solve()
    #print('\n'*30)
    print(('ILP lower bound is {}'.format(lp.solution.get_objective_value())))
    #for e in H.hyperedge_id_iterator():
        #if lp.solution.get_values(a(e)) == 1:
            #print(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e))

    #print('\n'*30)
    return lp

def build_initial_LP(lpfile,cuts,realvaluedcuts=False,heuristicsolution=[]):
    '''

    '''
    lp = cplex.Cplex()
    lp.read(lpfile)
    I = get_addable_cuts(cuts,realvaluedcuts=realvaluedcuts)
    for j in I:
        lp.linear_constraints.add(lin_expr = [j], senses = ['G'], rhs = [1], names = ['starter{0}'.format(j)])
    #print ('______________numconstraints____________')
    #print(len(I))
    if heuristicsolution != []:
        #The heuristic solution should be a list of edges
        ones = [a(c) for c in heuristicsolution]
        val = [1.0]*len(ones)
        eq = cplex.SparsePair(ind=ones,val=val)
        lp.MIP_starts.add(eq,lp.MIP_starts.effort_level.auto)
    return lp

def get_LP_lower_bound(lp):
    '''

    '''
    lp.solve()
    print(('LP had {} variables'.format(lp.variables.get_num())))
    print(('LP had {} non-binary variables'.format(len([a for a in lp.solution.get_values() if a > 0 and a < 1]))))
    #print('\n'*30)
    print(('LP lower bound is {}'.format(lp.solution.get_objective_value())))
    #print('\n'*30)

def get_LP_trim_order(edgelist,lp,edgedict):
    '''
    initially sets all edges to not be forced into the solution, then one by one forces them into the solution and sets their order as the LP value with them forced in the solution
    '''
    for e in edgelist:
        lp.variables.set_lower_bounds(a(e),0.0)
        
    for e in edgelist:
        lp.variables.set_lower_bounds(a(e),1.0)
        #lp.variables.set_upper_bounds(a(e),1.0)
        lp.solve()
        #write_and_print_vars(lp,H,nodeset,numsolsfound,verbose,outprefix)
        edgedict[e] = lp.solution.get_objective_value()
        lp.variables.set_lower_bounds(a(e),0.0)
    return edgedict

def get_fixed_stoichiometries(H):
    stoichiometries = {}
    for e in H.hyperedge_id_iterator():
        stoichiometries[e] = {}
        for v in H.get_hyperedge_tail(e):
            #These should be set to the max of outdeg - 1 and 1
            stoichiometries[e][v] = 1.0 / max(1,len(H.get_forward_star(v)))
        for v in H.get_hyperedge_head(e):
            stoichiometries[e][v] = 1.0
    return stoichiometries


def write_covering_constraints(H,ilpfile,sources):

    for v in H.node_iterator():
        if v not in sources:
            inedges = H.get_backward_star(v)
            outedges = H.get_forward_star(v)
            for e in outedges:
                ilpfile.write(' - {} '.format(a(e)))
                for f in inedges:
                    ilpfile.write('+ {} '.format(a(f)))

                ilpfile.write('>= 0\n')

def write_hyperpath_target_constraint(H,ilpfile,target):
    for e in H.get_backward_star(target):
        ilpfile.write(' + {}'.format(a(e)))
    ilpfile.write(' >= 1\n')

def write_nonuseless_edge_constraints(H,ilpfile,target):
    '''
    Need to not consider any hyperedge in the backward_star of the target
    '''
    tinedges = H.get_backward_star(target)
    for e in H.hyperedge_id_iterator():
        if e not in tinedges:
            outedges = set()
            for v in H.get_hyperedge_head(e):
                outedges.update(H.get_forward_star(v))
            for f in outedges:
                ilpfile.write(' + {} '.format(a(f)))
            ilpfile.write(' - {} >= 0\n'.format(a(e)))

def write_shortest_hyperpath_hybrid_ilp(H,ilpfilename,target,source,notchh=False,notc=False,nohh=False,notarget=False):
    ilpfile = open(ilpfilename,'w')

    ## write objective
    writeObjective(H,ilpfile,minimize=True)
    ilpfile.write('\nSubject To\n')
    if notarget == False:
        write_hyperpath_target_constraint(H,ilpfile,target)
    if notchh == False and notc == False:
        write_covering_constraints(H,ilpfile,[source])
    if notchh == False and nohh == False:
        write_nonuseless_edge_constraints(H,ilpfile,target)
    ilpfile.write('Binary\n')
    for e in H.hyperedge_id_iterator():
        ilpfile.write('{}\n'.format(a(e)))

    ilpfile.write('End\n')
    ilpfile.close()

def find_shortest_hyperpath_parallel(H,nodeset,ilpfilename,target,source,node_dict,headcuts=[],ilpchar='a',vertexheadcuts=[],nodb=False,notchh=False,noaugment=False,notc=False,nohh=False,notarget=False,returnP=False):
    '''
    Creates the ilp object for shortest hyperpath that includes factory flux constraints to make the target in quantity epsilon = 1, the new covering constraints and the new nonuseless edge constraints
    Inputs: H, hypergraph; nodeset, nodes in H; ilpfilename, the name of the file where the ilp will be written then subsequently read; target, the target; source, the source; headcuts, the headcuts that were calculated from the heuristic hyperpath
    Outputs: None
    Prints: The number of cutting plane iterations necessary to return a hyperpath, and the objective value for that final hyperpath
    Filewrites: For each iteration, the iteration number and the objective value should be output to a file.
    '''
    #TODO: Add the output solution to the output file, along with the flux values for each hyperedge

    write_shortest_hyperpath_hybrid_ilp(H,ilpfilename,target,source,notchh=notchh,notc=notc,nohh=nohh,notarget=notarget)
    
    #add in the head cuts
    if nodb == False:
        ilp = build_initial_LP(ilpfilename,headcuts)
    else:
        #This means exclude the head cuts becuase we only want the tail-covering and head-hitting cuts
        ilp = build_initial_LP(ilpfilename,[])
    #This is now an ilp object with the shortest hyperpath objective, the headcuts and factory constraints, and binary constraints on the hyperedge variables, so we just need to run cutting plane
    outprefix = 'reactomefachyperpathoutputs/hybridhyperpathilp{}.out'.format(ilpchar)
    if returnP == False:
        itrs,obj = solveILPcuttingplanesparallel(H,ilp,target,source,outprefix,nodeset,node_dict,numitrs = 10000000,ilpchar=ilpchar,headcuts=vertexheadcuts,noaugment=noaugment,nodb=nodb)
    else:
        itrs,obj,P = solveILPcuttingplanesparallel(H,ilp,target,source,outprefix,nodeset,node_dict,numitrs = 10000000,ilpchar=ilpchar,headcuts=vertexheadcuts,noaugment=noaugment,nodb=nodb,returnP=returnP)

    print('number of iterations: {}, objective value: {}'.format(itrs,obj))
    if returnP == True:
        return H,P

def find_shortest_hyperpath_with_new_inequalities(H,nodeset,ilpfilename,target,source,headcuts=[],ilpchar='a',heuristicpath=[],vertexheadcuts=[]):
    '''
    Creates the ilp object for shortest hyperpath that includes factory flux constraints to make the target in quantity epsilon = 1, the new covering constraints and the new nonuseless edge constraints
    Inputs: H, hypergraph; nodeset, nodes in H; ilpfilename, the name of the file where the ilp will be written then subsequently read; target, the target; source, the source; headcuts, the headcuts that were calculated from the heuristic hyperpath
    Outputs: None
    Prints: The number of cutting plane iterations necessary to return a hyperpath, and the objective value for that final hyperpath
    Filewrites: For each iteration, the iteration number and the objective value should be output to a file.
    '''
    #TODO: Add the output solution to the output file, along with the flux values for each hyperedge

    write_shortest_hyperpath_hybrid_ilp(H,ilpfilename,target,source)
    
    #add in the head cuts
    ilp = build_initial_LP(ilpfilename,headcuts,heuristicsolution=heuristicpath)
    #This is now an ilp object with the shortest hyperpath objective, the headcuts and factory constraints, and binary constraints on the hyperedge variables, so we just need to run cutting plane
    outprefix = 'reactomefachyperpathoutputs/hybridhyperpathilp{}.out'.format(ilpchar)
    itrs,obj = solveILPcuttingplanes(H,ilp,target,source,outprefix,nodeset,numitrs = 2500,ilpchar=ilpchar,headcuts=vertexheadcuts,heuristicobj=weight(H,heuristicpath))
    print('number of iterations: {}, objective value: {}'.format(itrs,obj))

def find_shortest_hyperpath_with_factory(H,nodeset,ilpfilename,target,source,headcuts=[],ilpchar='a',realvaluedcuts=False):
    '''
    Creates the ilp object for shortest hyperpath that includes factory flux constraints to make the target in quantity epsilon = 1, and keep conservation at each vertex.
    NOTE: we had to change the flux variables so they are all integers (otherwise the integer variables can be set by cut inequalities, but their flux variables are still zero).
    Inputs: H, hypergraph; nodeset, nodes in H; ilpfilename, the name of the file where the ilp will be written then subsequently read; target, the target; source, the source; headcuts, the headcuts that were calculated from the heuristic hyperpath
    Outputs: None
    Prints: The number of cutting plane iterations necessary to return a hyperpath, and the objective value for that final hyperpath
    Filewrites: For each iteration, the iteration number and the objective value should be output to a file.
    '''
    #TODO: Add the output solution to the output file, along with the flux values for each hyperedge

    write_shortest_hyperpath_factory_ilp(H,ilpfilename,target,source,realvaluedcuts=realvaluedcuts)
    
    #add in the head cuts, realvaluedcuts will make them over the flux variables instead of the binary variables
    ilp = build_initial_LP(ilpfilename,headcuts,realvaluedcuts=realvaluedcuts)
    #This is now an ilp object with the shortest hyperpath objective, the headcuts and factory constraints, and binary constraints on the hyperedge variables, so we just need to run cutting plane
    outprefix = 'reactomefachyperpathoutputs/fachyperpathilp{}.out'.format(ilpchar)
    itrs,obj = solveILPcuttingplanes(H,ilp,target,source,outprefix,nodeset,numitrs = 2500,ilpchar=ilpchar,realvaluedcuts=realvaluedcuts)
    print(itrs,obj)


def run_LP_heuristics(H,nodeset,lpfile,outprefix,numsols,target,source,headcuts,augmentedcuts,subopt=False,verbose=False,targetname='SUPERTARGET'):
    '''

    '''
    #print(cuts)
    edgedict = {}
    for e in H.hyperedge_id_iterator():
        edgedict[e] = 0
    begin = time.time()

    lp = build_initial_LP(lpfile,headcuts)
    lp2 = build_initial_LP(lpfile,augmentedcuts)
    lp3 = build_initial_LP(lpfile,[])
    logfile = open('lplog.log','w')
    lp.set_results_stream(logfile)
    lp.set_error_stream(logfile)
    lp2.set_results_stream(logfile)
    lp2.set_error_stream(logfile)
    lp3.set_results_stream(logfile)
    lp3.set_error_stream(logfile)
    #lp4.set_results_stream(logfile)
    #lp4.set_error_stream(logfile)
    #lp5.set_results_stream(logfile)
    #lp5.set_error_stream(logfile)
    #lp5.set_warning_stream(logfile)
    #lp6.set_error_stream(logfile)
    #lp6.set_warning_stream(logfile)
    S = [source]

    get_LP_lower_bound(lp)
    get_LP_lower_bound(lp2)
    lptime = time.time()
    print(('LP lower bound took {}'.format(lptime - begin)))

    fulledgedict = get_LP_trim_order([e for e in H.hyperedge_id_iterator()],lp2,edgedict)
    P,_ = trim(H,H.hyperedge_id_iterator(),target,source,fulledgedict)
    print(('trim of full edge list after forcing edges in LP gives {} with weight {}'.format(P,weight(H,P))))
    trimtime = time.time()
    print(('LP trim took {}'.format(trimtime - lptime)))

    #ilp = ILP_from_LP(lp2,H,nodeset)
    ilph = ILP_from_LP(lp,H,nodeset)
    print('ilp with head cuts')
    ilpa = ILP_from_LP(lp2,H,nodeset)
    print('ilp with augmented cuts')
    ilpe = ILP_from_LP(lp3,H,nodeset)
    ilptime = time.time()
    print(('ILP lower bound took {}'.format(ilptime - trimtime)))

    ilph.set_warning_stream(logfile)
    ilpe.set_warning_stream(logfile)
    ilpa.set_warning_stream(logfile)
    #ilph.set_warning_stream(logfile)
    solveILPcuttingplanes(H,ilpa,target,source,outprefix,nodeset,numitrs = 2500,targetname=targetname,ilpchar='a')
    #solveILPcuttingplanes(H,ilph,target,source,outprefix,nodeset,numitrs = 2500,targetname=targetname,ilpchar='h')
    #solveILPcuttingplanes(H,ilpe,target,source,outprefix,nodeset,numitrs = 2500,targetname=targetname,ilpchar='e')
    #numiterations = solveILPcuttingplanes(H,ilpb,target,source,outprefix,nodeset,numitrs = 25000,targetname=targetname,ilpchar='b')
    #fale = open('{}numiterations.txt'.format(outprefix),'w')
    #fale.write(str(numiterations))
    #fale.close()
    #print(numiterations)
    #cuttingplanestime = time.time()
    #print(('cutting planes took {}'.format(cuttingplanestime - ilptime)))

    LP_relaxation_heuristic(H,target,source,lp2)
    #relaxationtime = time.time()
    #print(('LP relaxation heuristic took {}'.format(relaxationtime - cuttingplanestime)))
    #print(('LP relaxation heuristic took {}'.format(relaxationtime - ilptime)))


def LP_relaxation_heuristic(H,target,source,lp):
    '''

    '''

    edgeremoved = set()
    S = set()
    scannablelist = set(H.hyperedge_id_iterator())
    LPedges = []
    while t_in_S(S,target) == False: 
        print('len path so far',len(edgeremoved))
        edgedistance = {}
        lp.solve()
        lowestobj = lp.solution.get_objective_value() #get the baseline objective
        for e in scannablelist:
            lp.variables.set_upper_bounds(a(e),0.0) #remove edge from LP
            try:
                lp.solve()
                diff = lp.solution.get_objective_value()
            except:
                #Means that the LP had no solution with this edge removed, so diff = infinity
                diff = 99999999999
            edgedistance[e] = diff
            lp.variables.set_upper_bounds(a(e),1.0) #put edge back in the LP
        biggestdiff = 0
        biggestdiffedge = None
        for e in edgedistance:
            if edgedistance[e] >= biggestdiff:
                biggestdiff = edgedistance[e]
                biggestdiffedge = e
        edgeremoved.add(biggestdiffedge)
        scannablelist.remove(biggestdiffedge)
        edgedistance.pop(e)
        ones = [a(edge) for edge in edgeremoved]
        LPedges = list(edgeremoved)
        S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
        lp.variables.set_lower_bounds(a(biggestdiffedge),1.0) #force it into the solution

    #Print no solution if the scannable set is empty and the sink is still not reachable from the source

    lpedgedict = {}
    for e in LPedges:
        lpedgedict[e] = 0
    lpedgedict = get_LP_trim_order(LPedges,lp,lpedgedict)
    P,_ = trim(H,LPedges,target,source,lpedgedict)
    #print('\n'*30)
    print(('LP relaxation heuristic gives the follwing path: {} with weight {}'.format(P,weight(H,P))))
    #print('\n'*30)

def print_path(H,node_dict,P):
    for e in P:
        print(e,[node_dict[v] for v in H.get_hyperedge_tail(e)],[node_dict[w] for w in H.get_hyperedge_head(e)])

def solveILPcuttingplanesparallel(H,ilp,target,source,outprefix,nodeset,node_dict,numitrs=100000,verbose=False,targetname='SUPERTARGET',ilpchar='2',headcuts=[],heuristicobj=None,noaugment=False,nodb=False,returnP=False):
    '''

    '''
    val1 = Value('d')
    val2 = Array('c',7000000)
    val3 = Array('c',700000)
    return_dict = [val1,val2,val3]
    p = Process(target=find_cuts_run_heuristic,args=(H,source,target,node_dict,return_dict,))
    p.start()
    numsolsfound = 1
    begin = time.time()
    iterationtime = begin
    numoptobjective = 0
    maxobj = None
    allvars = []
    S = [source]
    numiterations = 0
    ofile = open('{}ilpdata{}.txt'.format(targetname[-6:],ilpchar),'w')
    unaones = None
    headcutsadded = False
    #unaones = pkl.load(open('unaones.pkl','rb'))

    while t_in_S(S,target) == False and numiterations < numitrs:
        ## Solve ILP
        if not p.is_alive() and headcutsadded == False:
            headcutsadded = True
            #print('palive',p.is_alive(),return_dict[1][:])
            #We now unserialize the vertex head cuts
            #serialedgeheadcuts = []
            serialvertexheadcuts = []
            for ii in range(len(return_dict[1])):
                if return_dict[1][ii] != b'\x00':
                    serialvertexheadcuts.append(return_dict[1][ii].decode())
                else:
                    print('broke at {}'.format(ii))
                    break
            headcuts = unserialize_vertex_headcuts(serialvertexheadcuts)
            #print(edgeheadcuts)
            #headcuts = edge_to_vertex_cuts(H,edgeheadcuts)
            print('len headcuts after heuristic is: {}'.format(len(headcuts)))
            heuristicsolution = []
            for ii in range(len(return_dict[2])):
                if return_dict[1][ii] != b'\x00':
                    heuristicsolution.append(return_dict[2][ii].decode())
                else:
                    print('broke at {}'.format(ii))
                    break
            heuristicsolution = unserialize_path(heuristicsolution)
            print(heuristicsolution)
            heuristicobj = return_dict[0].value
            print(heuristicobj)
            hones = [a(c) for c in heuristicsolution]
            hval = [1.0]*len(hones)
            print(hones,hval)
            heq = cplex.SparsePair(ind=hones,val=hval)
            ilp.MIP_starts.add(heq,ilp.MIP_starts.effort_level.auto)
            headcutstuples = set()
            for cut in headcuts:
                headcutstuples.add(tuple(cut))
            edgeheadcuts = convert_vertex_cuts_to_edge_cuts(H,headcutstuples,node_dict)
            I = get_addable_cuts(edgeheadcuts)
            print('len addable cuts: {}'.format(len(I)))
            if nodb == False:
                for jj in I:
                    ilp.linear_constraints.add(lin_expr = [jj], senses = ['G'], rhs = [1], names = ['starter{0}'.format(jj)])
        elif p.is_alive():
            print('heuristic still running')
        if headcuts != [] and unaones != None and noaugment != True:
            sinkcutcount = 0
            sourcecutcount = 0
            #print(len(headcuts),'headcutlen')
            for hcount,h in enumerate(headcuts):
                #print(h)
                sinkcut = sinkcut_from_headcut(H,unaones,h)
                if source in sinkcut and targetname not in sinkcut:
                    sinkcutcount += 1
                    crossedges = get_new_constraint(sinkcut,H)
                    #print('got crossedges')
                    ilpcrossedges = [a(c) for c in crossedges]
                    val = [1.0]*len(ilpcrossedges)
                    eq = cplex.SparsePair(ind=ilpcrossedges,val=val)
                    ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%dheadcut%dsink' % (numiterations,hcount)])


                sourcecut = sourcecut_from_headcut(H,unaones,h)
                #print(sourcecut)
                if source in sourcecut and targetname not in sourcecut:
                    sourcecutcount += 1
                    crossedges = get_new_constraint(sourcecut,H)
                    #print('got crossedges')
                    ilpcrossedges = [a(c) for c in crossedges]
                    val = [1.0]*len(ilpcrossedges)
                    eq = cplex.SparsePair(ind=ilpcrossedges,val=val)
                    ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%dheadcut%dsource' % (numiterations,hcount)])
            print('cut counts:', sinkcutcount,sourcecutcount)

        else:
            crossedges = get_new_constraint(S,H)
            ones = [a(c) for c in crossedges]
            val = [1.0]*len(ones)
            eq = cplex.SparsePair(ind=ones,val=val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%d' % (numiterations)])

        ilp.solve()
        
        if ilp.solution.pool.get_num()>0:
            #print ('Solution Found.')
            #extract the edges from the solution
            objective = ilp.solution.pool.get_objective_value(0)
            print('objective',objective)
            ofile.write('iteration {} solval {}\n'.format(numiterations,objective))
            aedges = [a(e) for e in H.hyperedge_id_iterator()]
            variables = ilp.solution.pool.get_values(0,aedges)
            #print('variables',variables)
            #print('aedges',aedges)
            #ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
            #variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            #ones = [var for var in variables if variables[var] > .5 and var[0:3]=='a_e' and var.split('_')[1] not in nodeset]
            ones = [aedges[i] for i in range(len(aedges)) if variables[i] > .5]
            onesvals = [variables[i] for i in range(len(aedges)) if variables[i] > .5]
            ones.sort()
            #onesvals = [variables[one] for one in ones]
            #fones = [var for var in variables if variables[var] > .01 and var[0:3]=='f_e' and var.split('_')[1] not in nodeset]
            #fones.sort()
            #fonesvals = [variables[fone] for fone in fones]
            unaones = [una(one) for one in ones]
            #pkl.dump(unaones,open('unaones.pkl','wb'))
            #for u in unaones:
                #if u in crossedges:
                    #print(u,'was in crossedges')
                    #break
            #else:
                #print('no edge in the solution was in the prior crossedges')
            #if verify_set_ratios_factory(H,unaones) == False:
                #print('THIS IS NOT A FACTORY\n\n\n\n\n\n\n\n\n')
            print (ones,'ones')
            print(onesvals, 'onesvals')
            iterationtimeold = iterationtime
            iterationtime = time.time()
            print('time for iteration', iterationtime - iterationtimeold)
            #print(fones,'fones')
            #print(fonesvals, 'fonesvals')
            #for f in ones:
                #f = una(f)
                #print('\nedge',f)
                #for v in H.get_hyperedge_tail(f):
                    #print(v[-5:])
                #print('head')
                #if len(H.get_hyperedge_head(f)) < 300:
                    #for v in H.get_hyperedge_head(f):
                        #print(v[-5:])
                #else:
                    #print('lots')
            S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
            #print ('S',S,temp2)
            S.add(source)
        else:
            print ('Infeasible Solution. quitting.')
            break
        numsolsfound+=1
        numiterations+=1
        if heuristicobj != None and objective == heuristicobj:
            print('objective was equal to heuristic length, heuristic is optimal')
            if returnP == True:
                return numiterations, objective, unaones
            else:
                return numiterations, objective

    allvars = variables
    if p.is_alive():
        p.terminate()
        print('heuristic did not finish')
    print_path(H,node_dict,unaones)
    #print ('-'*20 + 'Cplex Output End' + '-'*20 + '\n')
    #print ('%d solutions found' % (numsolsfound-1))
    if returnP == True:
        return numiterations, objective, unaones
    else:
        return numiterations, objective

def solveILPcuttingplanes(H,ilp,target,source,outprefix,nodeset,numitrs=100000,verbose=False,targetname='SUPERTARGET',ilpchar='2',realvaluedcuts=False,headcuts=[],heuristicobj=None):
    '''

    '''
    numsolsfound = 1
    begin = time.time()
    iterationtime = begin
    numoptobjective = 0
    maxobj = None
    allvars = []
    S = [source]
    numiterations = 0
    ofile = open('{}ilpdata{}.txt'.format(targetname[-6:],ilpchar),'w')
    unaones = None
    #unaones = pkl.load(open('unaones.pkl','rb'))

    while t_in_S(S,target) == False and numiterations < numitrs:
        ## Solve ILP
        #NOTE: this gets what is now called a source constraint
        if headcuts != [] and unaones != None:
            sinkcutcount = 0
            sourcecutcount = 0
            for hcount,h in enumerate(headcuts):
                sinkcut = sinkcut_from_headcut(H,unaones,h)
                if source in sinkcut and targetname not in sinkcut:
                    sinkcutcount += 1
                    crossedges = get_new_constraint(sinkcut,H)
                    if realvaluedcuts == True:
                        ilpcrossedges = [fedge(c) for c in crossedges]
                    else:
                        ilpcrossedges = [a(c) for c in crossedges]
                    val = [1.0]*len(ilpcrossedges)
                    eq = cplex.SparsePair(ind=ilpcrossedges,val=val)
                    ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%dheadcut%dsink' % (numiterations,hcount)])


                sourcecut = sourcecut_from_headcut(H,unaones,h)
                if source in sourcecut and targetname not in sourcecut:
                    sourcecutcount += 1
                    crossedges = get_new_constraint(sourcecut,H)
                    if realvaluedcuts == True:
                        ilpcrossedges = [fedge(c) for c in crossedges]
                    else:
                        ilpcrossedges = [a(c) for c in crossedges]
                    val = [1.0]*len(ilpcrossedges)
                    eq = cplex.SparsePair(ind=ilpcrossedges,val=val)
                    ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%dheadcut%dsource' % (numiterations,hcount)])
            print('cut counts:', sinkcutcount,sourcecutcount)

        else:
            crossedges = get_new_constraint(S,H)
            if realvaluedcuts == True:
                ones = [fedge(c) for c in crossedges]
            else:
                ones = [a(c) for c in crossedges]
            val = [1.0]*len(ones)
            eq = cplex.SparsePair(ind=ones,val=val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%d' % (numiterations)])

        ilp.solve()
        
        if ilp.solution.pool.get_num()>0:
            #print ('Solution Found.')
            #extract the edges from the solution
            objective = ilp.solution.pool.get_objective_value(0)
            print('objective',objective)
            ofile.write('iteration {} solval {}\n'.format(numiterations,objective))
            aedges = [a(e) for e in H.hyperedge_id_iterator()]
            variables = ilp.solution.pool.get_values(0,aedges)
            #print('variables',variables)
            #print('aedges',aedges)
            #ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
            #variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            #ones = [var for var in variables if variables[var] > .5 and var[0:3]=='a_e' and var.split('_')[1] not in nodeset]
            ones = [aedges[i] for i in range(len(aedges)) if variables[i] > .5]
            onesvals = [variables[i] for i in range(len(aedges)) if variables[i] > .5]
            ones.sort()
            #onesvals = [variables[one] for one in ones]
            #fones = [var for var in variables if variables[var] > .01 and var[0:3]=='f_e' and var.split('_')[1] not in nodeset]
            #fones.sort()
            #fonesvals = [variables[fone] for fone in fones]
            unaones = [una(one) for one in ones]
            #pkl.dump(unaones,open('unaones.pkl','wb'))
            #for u in unaones:
                #if u in crossedges:
                    #print(u,'was in crossedges')
                    #break
            #else:
                #print('no edge in the solution was in the prior crossedges')
            #if verify_set_ratios_factory(H,unaones) == False:
                #print('THIS IS NOT A FACTORY\n\n\n\n\n\n\n\n\n')
            print (ones,'ones')
            print(onesvals, 'onesvals')
            iterationtimeold = iterationtime
            iterationtime = time.time()
            print('time for iteration', iterationtime - iterationtimeold)
            #print(fones,'fones')
            #print(fonesvals, 'fonesvals')
            #for f in ones:
                #f = una(f)
                #print('\nedge',f)
                #for v in H.get_hyperedge_tail(f):
                    #print(v[-5:])
                #print('head')
                #if len(H.get_hyperedge_head(f)) < 300:
                    #for v in H.get_hyperedge_head(f):
                        #print(v[-5:])
                #else:
                    #print('lots')
            S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
            #print ('S',S,temp2)
            S.add(source)
        else:
            print ('Infeasible Solution. quitting.')
            break
        numsolsfound+=1
        numiterations+=1
        if heuristicobj != None and objective == heuristicobj:
            print('objective was equal to heuristic length, heuristic is optimal')
            return numiterations, objective

    allvars = variables
    #print ('-'*20 + 'Cplex Output End' + '-'*20 + '\n')
    #print ('%d solutions found' % (numsolsfound-1))
    return numiterations, objective


def verify_set_ratios_factory(H,edges):
    intermediatemetabolites = set()
    for e in edges:
        intermediatemetabolites.update(H.get_hyperedge_tail(e))
        intermediatemetabolites.update(H.get_hyperedge_head(e))
    if 'SUPERTARGET' not in intermediatemetabolites:
        #This implies the target was not created
        return False
    for v in intermediatemetabolites:
        if 'SUPER' not in v:
            #we know fluxes were all set to 1, we can get what the stoichiometries are set to by 1 in head, 1/outdegree for tail
            incount = 0
            for e in H.get_backward_star(v):
                if e in edges:
                    incount += 1
            outcount = 0
            for e in H.get_forward_star(v):
                if e in edges:
                    outcount += 1
            inratio = 1
            outratio = 1.0 / max(1,len(H.get_forward_star(v)))

            if inratio * incount < outratio * outcount:
                print(v,incount,outcount, inratio, outratio)
                return False
    return True

################################

def solveILP(H,nodeset,lpfile,outprefix,numsols,subopt=False,verbose=False,printobjective=False):
    '''
    stuff
    '''
    print('\nSolving ILP...')

    print('\n' + '-'*20 + 'Cplex Output Start' + '-'*20)
    ilp = cplex.Cplex()
    ilp.parameters.threads.set(2)
    ilp.read(lpfile)

    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []
    objective = -1
    while(numsolsfound < numsols+1):
        ## Solve ILP
        print('-'*10 + 'Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n')
        ilp.solve()
        print('-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n')

        if ilp.solution.pool.get_num()>0:
            print('Solution Found.')

            objective = ilp.solution.pool.get_objective_value(0)
            if numsolsfound == 1:
                maxobj = objective
                numoptobjective+=1
                print('Max Objective of %d' % (objective))
            elif objective != maxobj and not subopt:
                print('Solution (obj=%d) does not have max objective of %d: quitting.' % (objective,maxobj))
                break
            elif objective != maxobj and subopt:
                print('Solution (obj=%d) does not have max objective of %d.' % (objective,maxobj))
            else:
                print('ANOTHER OPTIMAL OBJECTIVE=%d' % (objective))
                numoptobjective+=1

            # Incumbent solution is 0 in the pool:
            ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)

            # Get variables
            variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            allvars.append(variables)

            # Add constraint for this solution.
            # See http://www-01.ibm.com/support/docview.wss?uid=swg21399929
            #ones = [var for var in variables if variables[var] == 1]
            #zeros = [var for var in variables if variables[var] == 0]
            #ind = ones + zeros
            #val = [-1.0]*len(ones)+[1.0]*len(zeros)
            #eq = cplex.SparsePair(ind,val)
            #ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1-len(ones)], names = ['solution%d' % (numsolsfound)])

            # all we need is the ones.
            ones = [var for var in variables if variables[var] == 1 and var[0]=='a' and var.split('_')[1] not in nodeset]
            val = [1.0]*len(ones)
            eq = cplex.SparsePair(ones,val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(ones)-1], names = ['solution%d' % (numsolsfound)])
        else:
            print('Infeasible Solution. quitting.')
            break
        numsolsfound+=1

    print('-'*20 + 'Cplex Output End' + '-'*20 + '\n')
    print('%d solutions found' % (numsolsfound-1))
    if printobjective == True:
        print('Objective value is {0}'.format(objective))
        return numsolsfound-1,numoptobjective,allvars,objective
    else:
        return numsolsfound-1,numoptobjective,allvars

def solveILPoldcuttingplanes(H,nodeset,lpfile,outprefix,numsols,target,source,subopt=False,verbose=False):
    heuristicpath = heuristic_path(H,source,target)
    heuristiclist = {}
    for edge in H.hyperedge_id_iterator():
        if edge in heuristicpath:
            heuristiclist[a(edge)] = 1
        else:
            heuristiclist[a(edge)] = 0
        
    print(( 'the heuristic path:',heuristicpath))
    print( '\nSolving ILP...')

    print(( '\n' + '-'*20 + 'Cplex Output Start' + '-'*20))
    ilp = cplex.Cplex()
    ilp.read(lpfile)
    I = inequalities_from_path(H,heuristicpath,source,target)
    I = set(I)
    for j in I:
        #print 'herepringint' ,j)
        ilp.linear_constraints.add(lin_expr = [j], senses = ['G'], rhs = [1], names = ['starter{0}'.format(j)])

    return 5,5,heuristiclist
    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []
    S = [source]
    numiterations = 0
    while(t_in_S(S,target) == False):
        ## Solve ILP
        crossedges = get_new_constraint(S,H)
        print(( '_______crossedged:',crossedges))
        #addconstraint to the ILP
        ones = [a(c) for c in crossedges]
        val = [1.0]*len(ones)
        print(( len(crossedges),len(val)))
        eq = cplex.SparsePair(ind=ones,val=val)
        #still needs changing
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['iteration%d' % (numiterations)])
        print ('______________numconstraints____________')
        print((ilp.linear_constraints.get_num()))
        print(('-'*10, ilp.linear_constraints.get_rhs()))
        print(('-'*10 + 'Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        ilp.solve()
        print(('-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        
        if ilp.solution.pool.get_num()>0:
            print ('Solution Found.')
            #extract the edges from the solution
            objective = ilp.solution.pool.get_objective_value(0)
            ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
            variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            ones = [var for var in variables if variables[var] == 1 and var[0:3]=='a_e' and var.split('_')[1] not in nodeset]
            print((ones,'ones'))
            S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
            #print((S,temp2))
            S.add(source)
            
            #objective = ilp.solution.pool.get_objective_value(0)
            #if numsolsfound == 1:
                #maxobj = objective
                #numoptobjective+=1
                #print 'Max Objective of %d' % (objective))
            #elif objective != maxobj and not subopt:
                #print 'Solution (obj=%d) does not have max objective of %d: quitting.' % (objective,maxobj))
                #break
            #elif objective != maxobj and subopt:
                #print 'Solution (obj=%d) does not have max objective of %d.' % (objective,maxobj)))
            #else:
                #print 'ANOTHER OPTIMAL OBJECTIVE=%d' % (objective))
                #numoptobjective+=1

            # Incumbent solution is 0 in the pool:
            #ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)

            # Get variables 
            #variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            #allvars.append(variables)

            # Add constraint for this solution.
            # See http://www-01.ibm.com/support/docview.wss?uid=swg21399929
            #ones = [var for var in variables if variables[var] == 1]
            #zeros = [var for var in variables if variables[var] == 0]
            #ind = ones + zeros
            #val = [-1.0]*len(ones)+[1.0]*len(zeros)
            #eq = cplex.SparsePair(ind,val)
            #ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1-len(ones)], names = ['solution%d' % (numsolsfound)])

            # all we need is the ones.
            #ones = [var for var in variables if variables[var] == 1 and var[0]=='a' and var.split('_')[1] not in nodeset]
        else:
            print ('Infeasible Solution. quitting.')
            break
        numsolsfound+=1
        numiterations+=1

    allvars = variables
    print(('-'*20 + 'Cplex Output End' + '-'*20 + '\n'))
    print(('%d solutions found' % (numsolsfound-1)))
    return numsolsfound-1,numoptobjective,allvars

def solveILP_cheats(H,nodeset,lpfile,outprefix,numsols,cheat_set,subopt=False,verbose=False):
    print ('\nSolving ILP...')

    print(('\n' + '-'*20 + 'Cplex Output Start' + '-'*20))
    ilp = cplex.Cplex()
    ilp.read(lpfile)

    CHEAT_NUM = 5

    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []
    prevconstraint = None
    UPPER_BOUND = 100
    iteration = 0
    while(numsolsfound < numsols+1 and iteration < UPPER_BOUND):
        iteration+=1
        ## Solve ILP
        print(('-'*10 + 'ADDING CONSTRAINT FOR <%d CHEATS' % (CHEAT_NUM) + '-'*10 ))
        if prevconstraint:
            ilp.linear_constraints.delete(prevconstraint)

        # all we need is the ones.
        ones = ['a_'+c for c in cheat_set]
        val = [1.0]*len(ones)
        eq = cplex.SparsePair(ones,val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [CHEAT_NUM], names = ['cheat'])
        prevconstraint = ['cheat']
        CHEAT_NUM += 1

        print(('-'*10 + 'Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        ilp.solve()
        print(('-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        
        if ilp.solution.pool.get_num()>0:
            print ('Solution Found.')
            
            objective = ilp.solution.pool.get_objective_value(0)
            if numsolsfound == 1:
                maxobj = objective
                numoptobjective+=1
                print(('Max Objective of %d' % (objective)))
            elif objective != maxobj and not subopt:
                print(('Solution (obj=%d) does not have max objective of %d: quitting.' % (objective,maxobj)))
                break
            elif objective != maxobj and subopt:
                print(('Solution (obj=%d) does not have max objective of %d.' % (objective,maxobj)))
            else:
                print(('ANOTHER OPTIMAL OBJECTIVE=%d' % (objective)))
                numoptobjective+=1

            # Incumbent solution is 0 in the pool:
            ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)

            # Get variables 
            variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
            allvars.append(variables)

            # Add constraint for this solution.
            # See http://www-01.ibm.com/support/docview.wss?uid=swg21399929
            #ones = [var for var in variables if variables[var] == 1]
            #zeros = [var for var in variables if variables[var] == 0]
            #ind = ones + zeros
            #val = [-1.0]*len(ones)+[1.0]*len(zeros)
            #eq = cplex.SparsePair(ind,val)
            #ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1-len(ones)], names = ['solution%d' % (numsolsfound)])

            # all we need is the ones.
            ones = [var for var in variables if variables[var] == 1 and var[0]=='a' and var.split('_')[1] not in nodeset]
            val = [1.0]*len(ones)
            eq = cplex.SparsePair(ones,val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(ones)-1], names = ['solution%d' % (numsolsfound)])
            numsolsfound+=1

    if iteration == UPPER_BOUND:
        print ('Infeasible Solution. quitting.')
        
    print(('-'*20 + 'Cplex Output End' + '-'*20 + '\n'))
    print(('%d solutions found' % (numsolsfound-1)))
    return numsolsfound-1,numoptobjective,allvars

def solveILP_cutting_planes(H,nodeset,lpfile,outprefix,numsols,subopt=False,verbose=False,source='SUPERSOURCE',target='SUPERTARGET'):
    print ('\nSolving ILP...')

    print(('\n' + '-'*20 + 'Cplex Output Start' + '-'*20))
    ilp = cplex.Cplex()
    ilp.read(lpfile)

    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []

    #x = ilp.variables.get_names()
    x = ['a_%s' % n for n in H.get_hyperedge_id_set()]
    while(numsolsfound < numsols+1):
        ## Solve ILP
        print(('-'*10 + 'Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        done = False
        t = 0
        oldset = set()
        while not done:
            print(('....looking (iter %d)....' % (t)))
            ilp.solve()
            if ilp.solution.pool.get_num() == 0:
                done = True
            else:
                ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
                nodes,edges = visit_set(H,x,ilp,source)
                nodes.remove(source)
                print(('%d edges in solution: %s' % (len(edges),';'.join(edges))))
                print(('%d nodes in B-visit set: %s' % (len(nodes),';'.join(nodes))))
                if nodes == oldset:
                    sys.exit()
                oldset = nodes
                if target in nodes:    
                    done = True
                else: ## add new constraint.
                    spanning = [e for e in H.get_hyperedge_id_set() if \
                        len(nodes.intersection(H.get_hyperedge_tail(e))) > 0 and \
                        len(nodes.intersection(H.get_hyperedge_head(e))) < len(H.get_hyperedge_head(e))]
                    print(('%d spanning hyperedges' % (len(spanning))))
                    for e in spanning:
                        print(('  %s: %s --> %s' % (e,';'.join(H.get_hyperedge_tail(e)),';'.join(H.get_hyperedge_head(e)))))
                    ones = ['a_%s' % s for s in spanning]
                    val = [1.0]*len(ones)
                    eq = cplex.SparsePair(ones,val)
                    ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1], names = ['cutting_plane_constraint_%d' % (t)])
                    t+=1
        print(('-'*10 + 'Done Looking for Solution %d' % (numsolsfound) + '-'*10 + '\n'))
        
        if ilp.solution.pool.get_num() == 0:
                print ('Infeasible Solution. quitting.')
                break
        print ('Solution Found.')
            
        objective = ilp.solution.pool.get_objective_value(0)
        if numsolsfound == 1:
            maxobj = objective
            numoptobjective+=1
            print(('Max Objective of %d' % (objective)))
        elif objective != maxobj and not subopt:
            print(('Solution (obj=%d) does not have max objective of %d: quitting.' % (objective,maxobj)))
            break
        elif objective != maxobj and subopt:
            print(('Solution (obj=%d) does not have max objective of %d.' % (objective,maxobj)))
        else:
            print(('ANOTHER OPTIMAL OBJECTIVE=%d' % (objective)))
            numoptobjective+=1

        # Incumbent solution is 0 in the pool:
        ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)

        # Get variables 
        variables = getILPSolution(H,nodeset,outprefix,numsolsfound,objective,verbose)
        allvars.append(variables)

        # Add constraint for this solution.
        # See http://www-01.ibm.com/support/docview.wss?uid=swg21399929
        #ones = [var for var in variables if variables[var] == 1]
        #zeros = [var for var in variables if variables[var] == 0]
        #ind = ones + zeros
        #val = [-1.0]*len(ones)+[1.0]*len(zeros)
        #eq = cplex.SparsePair(ind,val)
        #ilp.linear_constraints.add(lin_expr = [eq], senses = ['G'], rhs = [1-len(ones)], names = ['solution%d' % (numsolsfound)])

        # all we need is the ones.
        ones = [var for var in variables if variables[var] == 1 and var[0]=='a' and var.split('_')[1] not in nodeset]
        val = [1.0]*len(ones)
        eq = cplex.SparsePair(ones,val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(ones)-1], names = ['solution%d' % (numsolsfound)])
            
        numsolsfound+=1

    print(('-'*20 + 'Cplex Output End' + '-'*20 + '\n'))
    print(('%d solutions found' % (numsolsfound-1)))
    return numsolsfound-1,numoptobjective,allvars

def visit_set(H,x,ilp,source):
    y = ilp.solution.get_values(x)
    hedges = [x[i].replace('a_','') for i in range(len(x)) if y[i]==1]
    HSUB = DirectedHypergraph()
    for hedge in hedges:
        HSUB.add_hyperedge(H.get_hyperedge_tail(hedge),H.get_hyperedge_head(hedge),weight=1)
    nodes,ignore,ignore = visit(HSUB,source)
    return nodes,hedges

################################
def getILPSolution(H,nodeset,outprefix,num,objective,verbose):
    #print ('\nGetting ILP Solution for Solution # %d in Pool' % (num))

    # parse xml
    xml = minidom.parse('%s-%d.sol' % (outprefix,num))
    cplexsol = xml.getElementsByTagName('CPLEXSolution')[0]
   
    # get variables
    elements = cplexsol.getElementsByTagName('variables')[0].getElementsByTagName('variable')
    variables = {}
    for v in elements:
        variables[str(v.getAttribute('name'))] = float(v.getAttribute('value'))

    out = open('%s-%d.variables' % (outprefix,num),'w')
    out.write('# Objective = %s\n' % (objective))
    out.write('#name\tval\trounded_val\n')

    if verbose == True:
        print ('VARIABLES:')
    
    numrounded = 0
    for v in variables:
        rounded = int(variables[v]+0.5)
        if verbose==True and variables[v] != 0:
            print((v,variables[v],rounded))
        out.write('%s\t%f\t%d\n' % (v,variables[v],rounded))
    out.close()
    
    #print (' wrote variables to file %s' % ('%s-%d.variables' % (outprefix,num)))

    return variables
    
def read_solution(solfile):

    alphas = {}
    pis = {}

    xmldoc = minidom.parse(solfile)
    varlist = xmldoc.getElementsByTagName('variable')
    for variable in varlist:
        name = variable.getAttribute('name').split('_')
        value = int(float(variable.getAttribute('value'))+0.5)

        if name[0] == 'a':
            alphas[name[1]] = value
        elif name[0] == 'pi':
            pis[name[1]] = value
        else:
            sys.exit('not tracking %s' % (name[0]))
    
    return alphas,pis

def a(name):
    return 'a_%s' % (name)

def o(name):
    return 'o_%s' % (name)

def pi(hnode,hedge):
    return 'pi_%s_%s' % (hnode,hedge)

def q(hedge):
    return 'q_%s' % (hedge)

def fedge(e):
    return('f_{}'.format(e))

def f(node1,node2):
    return 'f_%s_%s' % (node1,node2)

def una(name):
    return name[2:].strip()
