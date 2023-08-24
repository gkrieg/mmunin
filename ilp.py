import cplex
from multiprocessing import Process, Value, Array, Event

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit,visit

from pathheuristic import weight
from cutfinder import find_cuts_run_heuristic, unserialize_vertex_headcuts, unserialize_path, convert_vertex_cuts_to_edge_cuts

def writeObjective(H,out,minimize=True):
    if minimize:
        out.write('Minimize\n')
    else:
        out.write('Maximize\n')

    for hedge in H.hyperedge_id_iterator():
        out.write(' + %d %s' % (H.get_hyperedge_attribute(hedge,'weight'),a(hedge)))
    out.write('\n')
    return out

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
    nodeset = set()
    H2 = DirectedHypergraph()
    H2.add_nodes(H.get_node_set())
    for edge in ones:
        H2.add_hyperedge(H.get_hyperedge_tail(una(edge)),H.get_hyperedge_head(una(edge)))
    
    nodeset = b_visit(H2,source)
    #print ('nodeset size:',len(nodeset[0]), nodeset[0])
    return nodeset





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

            
def build_initial_LP(lpfile,cuts,realvaluedcuts=False,heuristicsolution=[]):
    '''

    '''
    lp = cplex.Cplex()
    lp.read(lpfile)
    I = get_addable_cuts(cuts,realvaluedcuts=realvaluedcuts)
    for j in I:
        lp.linear_constraints.add(lin_expr = [j], senses = ['G'], rhs = [1], names = ['starter{0}'.format(j)])
    if heuristicsolution != []:
        #The heuristic solution should be a list of edges
        ones = [a(c) for c in heuristicsolution]
        val = [1.0]*len(ones)
        eq = cplex.SparsePair(ind=ones,val=val)
        lp.MIP_starts.add(eq,lp.MIP_starts.effort_level.auto)
    return lp

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

def find_shortest_hyperpath_parallel(H,nodeset,ilpfilename,target,source,node_dict,headcuts=[],vertexheadcuts=[],nodb=False,notchh=False,noaugment=False,notc=False,nohh=False,notarget=False,returnP=False,outputfile=f'results/results.txt'):
    '''
    Creates the ilp object for shortest hyperpath that includes factory flux constraints to make the target in quantity epsilon = 1, the new covering constraints and the new nonuseless edge constraints
    Inputs: H, hypergraph; nodeset, nodes in H; ilpfilename, the name of the file where the ilp will be written then subsequently read; target, the target; source, the source; headcuts, the headcuts that were calculated from the heuristic hyperpath
    Outputs: None
    Prints: The number of cutting plane iterations necessary to return a hyperpath, and the objective value for that final hyperpath
    Filewrites: For each iteration, the iteration number and the objective value should be output to a file.
    '''

    write_shortest_hyperpath_hybrid_ilp(H,ilpfilename,target,source,notchh=notchh,notc=notc,nohh=nohh,notarget=notarget)
    
    #add in the head cuts
    if nodb == False:
        ilp = build_initial_LP(ilpfilename,headcuts)
    else:
        #This means exclude the head cuts becuase we only want the tail-covering and head-hitting cuts
        ilp = build_initial_LP(ilpfilename,[])
    #This is now an ilp object with the shortest hyperpath objective, the headcuts and factory constraints, and binary constraints on the hyperedge variables, so we just need to run cutting plane
    if returnP == False:
        itrs,obj = solveILPcuttingplanesparallel(H,ilp,target,source,nodeset,node_dict,numitrs = 10000000,headcuts=vertexheadcuts,outputfile=outputfile)
    else:
        itrs,obj,P = solveILPcuttingplanesparallel(H,ilp,target,source,nodeset,node_dict,numitrs = 10000000,headcuts=vertexheadcuts,returnP=returnP,outputfile=outputfile)

    print('number of iterations: {}, objective value: {}'.format(itrs,obj))
    if returnP == True:
        return H,P

def solveILPcuttingplanesparallel(H,ilp,target,source,nodeset,node_dict,numitrs=100000,verbose=False,targetname='SUPERTARGET',headcuts=[],heuristicobj=None,noaugment=False,nodb=False,returnP=False,outputfile=f'results/results.txt'):
    '''

    '''
    output = open(outputfile,'w')
    val1 = Value('d')
    val2 = Array('c',7000000)
    val3 = Array('c',700000)
    return_dict = [val1,val2,val3]
    event = Event()
    p = Process(target=find_cuts_run_heuristic,args=(H,source,target,node_dict,return_dict,event,))
    p.start()
    numsolsfound = 1
    numoptobjective = 0
    maxobj = None
    allvars = []
    S = [source]
    numiterations = 0
    unaones = None
    headcutsadded = False

    while t_in_S(S,target) == False and numiterations < numitrs:
        ## Solve ILP
        if not p.is_alive() and headcutsadded == False:
            headcutsadded = True
            #We now unserialize the vertex head cuts
            serialvertexheadcuts = []
            for ii in range(len(return_dict[1])):
                if return_dict[1][ii] != b'\x00':
                    serialvertexheadcuts.append(return_dict[1][ii].decode())
                else:
                    print('broke at {}'.format(ii))
                    break
            headcuts = unserialize_vertex_headcuts(serialvertexheadcuts)
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
            #extract the edges from the solution
            objective = ilp.solution.pool.get_objective_value(0)
            print('objective',objective)
            output.write(f'objective {objective}\n')
            aedges = [a(e) for e in H.hyperedge_id_iterator()]
            variables = ilp.solution.pool.get_values(0,aedges)
            ones = [aedges[i] for i in range(len(aedges)) if variables[i] > .5]
            onesvals = [variables[i] for i in range(len(aedges)) if variables[i] > .5]
            ones.sort()
            unaones = [una(one) for one in ones]
            print (ones,'ones')
            output.write(f'ones {ones}\n')
            print(onesvals, 'onesvals')
            output.write(f'onesvals {onesvals}\n')
            S,temp1,temp2,temp3 = reachability_from_edges(H,ones,source)
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
        event.set()
        p.join()
        print('heuristic did not finish')
    path_printed = print_path(H,node_dict,unaones)
    output.write(path_printed)
    if returnP == True:
        return numiterations, objective, unaones
    else:
        return numiterations, objective

def t_in_S(S,target):
    #here is some code about checking if the sink is in S
    #returns a boolean
    for node in S:
        if node == target:
            return True
    return False

def print_path(H,node_dict,P):
    ret_list = []
    for e in P:
        tail = [node_dict[v] for v in H.get_hyperedge_tail(e)]
        head = [node_dict[w] for w in H.get_hyperedge_head(e)]
        print(e,tail,head)
        ret_list.append(f"edge {e}, tail {tail}, head {head}")
    return '\n'.join(ret_list)

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
