#include <unistd.h> // getopt 
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>

#include <fstream> //file output
#include <sstream> //string operation

#include <stdio.h>

#include <thread>
#include <unistd.h> //sleep for secs
#include <time.h>

using namespace std;

//
// Channels
//
struct Cross_Paths {
   // Communication node pair index going through a channel
   vector<int> pair_index;
   // Communication flow index going through a channel
   vector<int> flow_index;   
   // list of ID in Pair for a channel
   vector<int> assigned_list;
   // list of h_dst in Pair for a channel
   vector<int> assigned_dst_list;
   // If all IDs are assigned
   bool Valid;
   // initialize
   Cross_Paths():Valid(false){}
   // routing table
   vector<int> routing_table;  // input port, slot number, src node, dst node, ...

   bool operator == (const Cross_Paths &a) const {
      return Valid == a.Valid && pair_index.size() == a.pair_index.size();
   }
   bool operator < (const Cross_Paths &a) const {
	 if (Valid == a.Valid)
	    return pair_index.size() < a.pair_index.size();
	 else if (Valid == false)
	    return false; 
	 else
	    return true;
   }
};

//
// Flow
//
struct Flow {   
   vector<int> pairs_id;
   // specified flow id
   int id; 
   // channels
   vector<int> channels;
   // assigned time slot #
   int ID;
};

//
// Communication node pair
//
struct Pair {
   // pair id (unique, from 0, 1, 2, ...)  
   int pair_id;
   // flow id (not unique)
   int flow_id; 
   // channels
   vector<int> channels;
   // source and destination
   int src; int dst;
   int h_src; int h_dst;
   // assigned ID (different from pair_id)
   int ID;
   // ID is valid or not
   bool Valid;
   // the distance between src and dst
   int hops;
   // initialize
   Pair(int s, int d, int h_s, int h_d): 
	src(s),dst(d),h_src(h_s),h_dst(h_d),flow_id(-1),ID(-1),Valid(false),hops(-1){}

   bool operator < (Pair &a){
	 if (Valid == a.Valid)
	    return hops < a.hops;
	 else if (Valid == false)
	    return false; 
	 else
	    return true;
   }
};

//
// Job
//
struct Job {   
   int time_submit; // submit time in workload
   int time_run; // runtime in workload
   int num_nodes; // number of nodes required
   //vector<int> src_dst_pair; // pairs of src and dst in workload
   //vector<int> src_dst_pair_m; // pairs of src and dst after dispatched
   vector<Pair> src_dst_pairs;
   vector<Pair> src_dst_pairs_m;
   int job_id; // job id (unique)
   time_t time_submit_r; // real submit time
   time_t time_dispatch_r; // real dispatch time
   vector<int> nodes; // occupied system nodes
   vector<Flow> flows; // one flow consists of one or multiple pairs
};

bool LSF(Job a, Job b) { return (a.num_nodes > b.num_nodes); } //large size first
bool SSF(Job a, Job b) { return (a.num_nodes < b.num_nodes); } //small size first
bool LRF(Job a, Job b) { return (a.time_run > b.time_run); } //large runtime first
bool SRF(Job a, Job b) { return (a.time_run < b.time_run); } //small runtime first

// ########################################## //
// ##  C program for Dijkstra's single source shortest path algorithm. 
// ##  The program is for adjacency matrix representation of the graph. ## //
// ##  Modified based on https://www.geeksforgeeks.org/printing-paths-dijkstras-shortest-path-algorithm/  ## //
// ########################################## //

// A utility function to find the vertex with minimum distance value, from the set of vertices not yet included in shortest path tree 
int minDistance(int dist[], bool sptSet[], int V) 
{    
    // Initialize min value 
    int min = 10000, min_index; 
  
    for (int v = 0; v < V; v++) 
        if (sptSet[v] == false && dist[v] <= min) {
                min = dist[v]; 
                min_index = v; 
        }
        
    return min_index; 
} 

// Function to print shortest path from source to j using parent array 
void printPath(int parent[], int j, int src, int dst, vector< vector<int> > &pair_path, int V) 
{ 
    // Base Case : If j is source 
    if (parent[j] == -1) 
        return; 
  
    printPath(parent, parent[j], src, dst, pair_path, V); 
  
    //printf("%d ", j); 
    pair_path[src*V+dst].push_back(j);
} 

// A utility function to print the constructed distance array 
int printSolution(int dist[], int V, int parent[], int src, vector< vector<int> > &pair_path) 
{ 
    //int src = 0; 
    //printf("Vertex\t Distance\tPath"); 
    for (int i = 0; i < V; i++) 
    { 
        if (i != src) {
                //printf("\n%d -> %d \t\t %d\t\t%d ", src, i, dist[i], src); 
                int dst = i;
                printPath(parent, i, src, dst, pair_path, V); 
        }
    } 
}

// Funtion that implements Dijkstra's single source shortest path algorithm for a graph represented using adjacency matrix representation 
void dijkstra(int V, vector<int> graph, int src, vector< vector<int> > &pair_path) 
{ 
      
    // The output array. dist[i] will hold the shortest distance from src to i 
    int dist[V];  
  
    // sptSet[i] will true if vertex i is included / in shortest path tree or shortest distance from src to i is finalized 
    bool sptSet[V]; 
  
    // Parent array to store shortest path tree 
    int parent[V]; 
  
    // Initialize all distances as INFINITE and stpSet[] as false 
    for (int i = 0; i < V; i++) 
    { 
        parent[i] = -1; 
        dist[i] = 10000; 
        sptSet[i] = false; 
    } 
  
    // Distance of source vertex from itself is always 0 
    dist[src] = 0; 
  
    // Find shortest path for all vertices 
    for (int count = 0; count < V - 1; count++) 
    { 
        // Pick the minimum distance vertex from the set of vertices not yet processed. u is always equal to src in first iteration. 
        int u = minDistance(dist, sptSet, V); 
  
        // Mark the picked vertex as processed 
        sptSet[u] = true; 
  
        // Update dist value of the adjacent vertices of the picked vertex. 
        for (int v = 0; v < V; v++) 
  
            // Update dist[v] only if is not in sptSet, there is an edge from u to v, and total weight of path from src to v through u is smaller than current value of dist[v] 
            if (!sptSet[v] && graph[u*V+v] && dist[u] + graph[u*V+v] < dist[v]) 
            { 
                parent[v] = u; 
                dist[v] = dist[u] + graph[u*V+v]; 
            }  
    } 
  
    // print the constructed distance array 
    printSolution(dist, V, parent, src, pair_path); 
}

void check_state(vector<int> & sys, vector<Job> & queue, bool & all_submitted, bool & all_finished){
    while (true) {
        int ava = 0;
        for (int i=0; i<sys.size(); i++){
            ava += sys[i];
        }
        if (ava == 0 && queue.size() == 0 && all_submitted == true){
            all_finished = true;
            break;
        }
    }
}

//job submit thread
void submit_jobs(vector<Job> & all_jobs, vector<Job> & queue, int queue_policy, bool & all_submitted){
        queue.push_back(all_jobs[0]);
        all_jobs[0].time_submit_r = time(NULL);
        int current_submit = all_jobs[0].time_submit;
        //cout << "Job " << all_jobs[0].job_id << " is submitted" << endl;
        printf("Job %d is submitted\n", all_jobs[0].job_id); 
        for (int i=1; i<all_jobs.size(); i++){
                int submit_interval = all_jobs[i].time_submit - current_submit;
                if (submit_interval > 0){
                        sleep(submit_interval);
                        queue.push_back(all_jobs[i]);
                        current_submit = all_jobs[i].time_submit;
                        if (queue.size()>2){
                                if (queue_policy == 1){ //bf
                                        sort(queue.begin()+1, queue.end(), LSF);
                                }
                                else if (queue_policy == 2){ //sf
                                        sort(queue.begin()+1, queue.end(), SSF);
                                } 
                                else if (queue_policy == 3){ //rlf
                                        sort(queue.begin()+1, queue.end(), LRF);
                                } 
                                else if (queue_policy == 4){ //rsf
                                        sort(queue.begin()+1, queue.end(), SRF);
                                } 
                        }  
                        all_jobs[i].time_submit_r = time(NULL);
                        //cout << "Job " << all_jobs[i].job_id << " is submitted" << endl;  
                        printf("Job %d is submitted\n", all_jobs[i].job_id);                                                                                       
                }
                else{
                        //cout << "Job " << all_jobs[i].job_id << " is not submitted due to errorous submit time" << endl;
                        printf("Job %d is not submitted due to errorous submit time\n", all_jobs[i].job_id);
                }
                if (i == all_jobs.size()-1){
                        all_submitted = true;
                        cout << "All jobs have been submitted" << endl;
                }
        }
}

void greeting(int a, vector<int> & b) {
    cout << "Hello multithread!" << a << endl;
    b.push_back(4);
}

vector<Pair> src_dst_pairs_map(vector<Pair> src_dst_pairs, vector<int> nodes, int topo, int hosts){
        vector<Pair> src_dst_pairs_m;
        vector<int> temp;
        for(int i=0; i<src_dst_pairs.size(); i++){
                //if (i%5 < 2) temp.push_back(src_dst_pair[i]); // i%5 == 0 --> src, i%5 == 1 --> dst, i%5 == 2 --> flow id, i%5 == 3 --> pair id, i%5 == 4 --> slot #
                temp.push_back(src_dst_pairs[i].h_src);
                temp.push_back(src_dst_pairs[i].h_dst);
        }
        sort(temp.begin(), temp.end());
        temp.erase(unique(temp.begin(), temp.end()), temp.end());
        for(int i=0; i<src_dst_pairs.size(); i++){
                // if (i%5 < 2){
                //         for(int j=0; j<temp.size(); j++){
                //                 if(src_dst_pair[i] == temp[j]){
                //                         src_dst_pair_m.push_back(nodes[j]);
                //                 }
                //         }
                // }
                // else{
                //         src_dst_pair_m.push_back(src_dst_pair[i]);
                // }
                int h_src_m = -1;
                int h_dst_m = -1;
                for(int j=0; j<temp.size(); j++){
                        if(src_dst_pairs[i].h_src == temp[j]){
                                h_src_m = nodes[j];
                        }
                        if(src_dst_pairs[i].h_dst == temp[j]){
                                h_dst_m = nodes[j];
                        }      
                        if (h_src_m != -1 and h_dst_m != -1){
                                int src_m = h_src_m/hosts;
                                int dst_m = h_dst_m/hosts;
                                if(topo == 2){ // fat-tree
                                        src_m = h_src_m;
                                        dst_m = h_dst_m;   
                                }
                                Pair pair(src_m,dst_m,h_src_m,h_dst_m);
                                pair.pair_id = src_dst_pairs[i].pair_id;
                                pair.flow_id = src_dst_pairs[i].flow_id;
                                src_dst_pairs_m.push_back(pair);  
                                break;                               
                        }                  
                }
        }
        return src_dst_pairs_m;
}

void update_after_release(Job job, vector<Cross_Paths> & Crossing_Paths)
{
    for (int p=0; p<job.src_dst_pairs_m.size(); p++){
       //update Crossing_Paths
       for (int i=0; i<job.src_dst_pairs_m[p].channels.size(); i++){
               vector<int> v = Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].pair_index;
               v.erase(remove(v.begin(), v.end(), job.src_dst_pairs_m[p].pair_id), v.end());
               v = Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].flow_index;
               v.erase(remove(v.begin(), v.end(), job.job_id*1000+job.src_dst_pairs_m[p].flow_id), v.end());  
               v = Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].assigned_list;
               v.erase(remove(v.begin(), v.end(), job.flows[job.job_id*1000+job.src_dst_pairs_m[p].flow_id].ID), v.end()); 
               v = Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].assigned_dst_list;
               v.erase(remove(v.begin(), v.end(), job.src_dst_pairs_m[p].h_dst), v.end());                                             
               for (int j=4; j<Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table.size(); j=j+5){
                       if (Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j] == job.src_dst_pairs_m[p].pair_id){
                               Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j-4] = -1;
                               Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j-3] = -1;
                               Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j-2] = -1;
                               Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j-1] = -1;
                               Crossing_Paths[job.src_dst_pairs_m[p].channels[i]].routing_table[j] = -1;
                       }
               }
       }
    }
}

void release_nodes(Job job, vector<int> & sys, vector<Cross_Paths> & Crossing_Paths)
{
       sleep(job.time_run);

       for (int n=0; n<job.nodes.size(); n++){
               sys[job.nodes[n]] = 0; //released
       }       

       update_after_release(job, Crossing_Paths);

       //cout << "Job " << job.job_id << " is finished" << endl;  
       printf("Job %d is finished\n", job.job_id); 
}

void route_job(int Topology, Job & job, int degree, int switch_num, int Host_Num, vector<Cross_Paths> & Crossing_Paths, int dimension, int array_size, int & hops, int & before_hops, int PORT, int node_num, 
                int groups, int group_switch_num, vector<int> & Switch_Topo, int inter_group, vector<int> topo_file, vector<int> topo_sws_uni, int & src, int & dst, int & h_src, int & h_dst, int & ct){
   if (Topology == 0 || Topology == 1) //mesh or torus
   {
        for (int i=0; i<job.src_dst_pairs_m.size(); i++){
            h_src = job.src_dst_pairs_m[i].h_src;
            h_dst = job.src_dst_pairs_m[i].h_dst;
            src = job.src_dst_pairs_m[i].src;
            dst = job.src_dst_pairs_m[i].dst;
            ct = job.src_dst_pairs_m[i].pair_id;

            bool wrap_around_x = false;
            bool wrap_around_y = false;
            bool wrap_around_z = false; //3D
            bool wrap_around_a = false; //4D

            //#######################//
            // e.g. 3D
            // (switch port: 0 (not used), 1 +x, 2 -x, 3 -y, 4 +y, 5 -z, 6 +z, 7 localhost->switch, 8 switch-> localhost)
            // switch port:1 +x, 2 -x, 3 +y, 4 -y, 5 -z, 6 +z,
                    // from 7 to (6+Host_Num) localhost->switch,
                    // from (7+Host_Num) to (6+Host_Num*2) switch-> localhost
            //#######################//

            // channel <-- node pair ID, node pair <-- channel ID 
            //Pair tmp_pair(src,dst,h_src,h_dst);  
            //pairs.push_back(tmp_pair);
            int t = src*(degree+1+2*Host_Num)+degree+1+h_src%Host_Num;
            // int t = (dst%2==1 && Topology==1) ? 
            //         Vch*src*(degree+1+2*Host_Num)+(degree+1+2*Host_Num)+degree+1+h_src%Host_Num
            //         : Vch*src*(degree+1+2*Host_Num)+degree+1+h_src%Host_Num;
            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
            job.src_dst_pairs_m[i].channels.push_back(t);  // node pair <-- channel ID  

            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

            //pairs[ct].pair_id = ct; 
            int delta_x, delta_y, delta_z, delta_a, current, src_xy, dst_xy, src_xyz, dst_xyz; // 2D, 3D, 4D
            if (dimension == 2 || dimension == 1){ //2D, 1D
                    src_xy = src; 
                    dst_xy = dst; 
            }
            if (dimension == 3){ //3D
                    src_xy = src%int(pow(array_size,2)); 
                    dst_xy = dst%int(pow(array_size,2)); 
            }
            if (dimension == 4){ //4D
                    src_xyz = src%int(pow(array_size,3)); 
                    dst_xyz = dst%int(pow(array_size,3)); 
                    src_xy = src_xyz%int(pow(array_size,2)); 
                    dst_xy = dst_xyz%int(pow(array_size,2));   
            }
            switch (Topology){
            case 0: //mesh
                    //  if (dimension == 2){ //2D
                    //         delta_x = dst%array_size - src%array_size;
                    //         delta_y = dst/array_size - src/array_size;
                    //  }
                    if (dimension == 3){ //3D
                            delta_z = dst/int(pow(array_size,2)) - src/int(pow(array_size,2)); 
                            // delta_x = dst_xy%array_size - src_xy%array_size; 
                            // delta_y = dst_xy/array_size - src_xy/array_size; 
                    }
                    if (dimension == 4){ //4D
                    delta_a = dst/int(pow(array_size,3)) - src/int(pow(array_size,3));
                    delta_z = dst_xyz/int(pow(array_size,2)) - src_xyz/int(pow(array_size,2)); 
                    //        delta_x = dst_xy%array_size - src_xy%array_size; 
                    //        delta_y = dst_xy/array_size - src_xy/array_size; 
                    }
                    //4D, 3D, 2D, 1D(delta_y=0)
                    delta_x = dst_xy%array_size - src_xy%array_size; 
                    delta_y = dst_xy/array_size - src_xy/array_size; 
                    current = src; 
                    break;

            case 1: // torus
                    if (dimension == 4){ //4D
                            delta_a = dst/int(pow(array_size,3)) - src/int(pow(array_size,3));
                            if ( delta_a < 0 && abs(delta_a) > array_size/2 ) {
                            //delta_a = -( delta_a + array_size/2);
                            delta_a = delta_a + array_size;
                                    wrap_around_a = true;		
                            } else if ( delta_a > 0 && abs(delta_a) > array_size/2 ) {
                            //delta_a = -( delta_a - array_size/2);
                            delta_a = delta_a - array_size;
                                    wrap_around_a = true;		
                            }
                            delta_z = dst_xyz/int(pow(array_size,2)) - src_xyz/int(pow(array_size,2));
                            if ( delta_z < 0 && abs(delta_z) > array_size/2 ) {
                            //delta_z = -( delta_z + array_size/2);
                            delta_z = delta_z + array_size;
                                    wrap_around_z = true;		
                            } else if ( delta_z > 0 && abs(delta_z) > array_size/2 ) {
                            //delta_z = -( delta_z - array_size/2);
                            delta_z = delta_z - array_size;
                                    wrap_around_z = true;		
                            }
                    }
                    if (dimension == 3){ //3D
                            delta_z = dst/int(pow(array_size,2)) - src/int(pow(array_size,2));
                            if ( delta_z < 0 && abs(delta_z) > array_size/2 ) {
                            //delta_z = -( delta_z + array_size/2);
                            delta_z = delta_z + array_size;
                                    wrap_around_z = true;		
                            } else if ( delta_z > 0 && abs(delta_z) > array_size/2 ) {
                            //delta_z = -( delta_z - array_size/2);
                            delta_z = delta_z - array_size;
                                    wrap_around_z = true;		
                            }
                    }
                    //4D, 3D, 2D, 1D(delta_y=0)
                    delta_x = dst_xy%array_size - src_xy%array_size;
                    if ( delta_x < 0 && abs(delta_x) > array_size/2 ) {
                    //delta_x = -( delta_x + array_size/2);
                    delta_x = delta_x + array_size;
                            wrap_around_x = true;		
                    } else if ( delta_x > 0 && abs(delta_x) > array_size/2 ) {
                    //delta_x = -( delta_x - array_size/2);
                    delta_x = delta_x - array_size;
                            wrap_around_x = true;		
                    }
                    delta_y = dst_xy/array_size - src_xy/array_size;
                    if ( delta_y < 0 && abs(delta_y) > array_size/2 ) {
                    //delta_y = -( delta_y + array_size/2);
                    delta_y = delta_y + array_size;
                            wrap_around_y = true;		
                    } else if ( delta_y > 0 && abs(delta_y) > array_size/2 ) {
                    //delta_y = -( delta_y - array_size/2);
                    delta_y = delta_y - array_size;
                            wrap_around_y = true;		
                    }
                    current = src; 
                    break;
            default:
                    cerr << "Please select -t0, or -t1 option" << endl;
                    exit(1);
                    break;
            }

            if (dimension == 4) //4D
            {
                    job.src_dst_pairs_m[i].hops = abs(delta_x) + abs(delta_y)+ abs(delta_z) + abs(delta_a);
            }
            if (dimension == 3) //3D
            {
                    job.src_dst_pairs_m[i].hops = abs(delta_x) + abs(delta_y)+ abs(delta_z);
            }
            if (dimension == 2 || dimension == 1) //2D, 1D(delta_y=0)
            {
                    job.src_dst_pairs_m[i].hops = abs(delta_x) + abs(delta_y);
            }

            if (dimension == 4){ //4D routing
                    if (delta_a > 0){
                            while ( delta_a != 0 ){  //-a
                            int t = current * (degree+1+2*Host_Num) + 7;
                            // int t = (wrap_around_a) ? Vch*current*(degree+1+2*Host_Num)+7+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 7;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == array_size-1) { 
                            if ( current >= array_size*array_size*array_size*(array_size-1)) {
                            wrap_around_a = false;
                            current = current - (array_size -1)*array_size*array_size*array_size;
                            } else current += array_size*array_size*array_size; 
                            delta_a--;
                            hops++;
                            }
                    } else if (delta_a < 0){
                            while ( delta_a != 0 ){  //+a
                            int t = current * (degree+1+2*Host_Num) + 8;
                            // int t = (wrap_around_a) ? Vch*current*(degree+1+2*Host_Num)+8+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 8;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == 0 ) { 
                            if ( current < array_size*array_size*array_size) {
                            wrap_around_a = false;
                            current = current + (array_size -1)*array_size*array_size*array_size;
                            } else current -= array_size*array_size*array_size;
                            hops++;
                            delta_a++;
                            }
                    }
                    
                    if (delta_a != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    } 

                    if (delta_z > 0){
                            while ( delta_z != 0 ){ // -z
                            int t = current * (degree+1+2*Host_Num) + 5;
                            // int t = (wrap_around_z) ? Vch*current*(degree+1+2*Host_Num)+5+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 5;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == array_size-1) { 
                            if ( (current%(array_size*array_size*array_size)) >= array_size*array_size*(array_size-1)) { 
                            wrap_around_z = false;
                            current = current - (array_size -1)*array_size*array_size;
                            } else current += array_size*array_size; 
                            delta_z--;
                            hops++;
                            }
                    } else if (delta_z < 0){
                            while ( delta_z != 0 ){ // +z
                            int t = current * (degree+1+2*Host_Num) + 6;
                            // int t = (wrap_around_z) ? Vch*current*(degree+1+2*Host_Num)+6+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 6;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == 0 ) { 
                            if ( (current%(array_size*array_size*array_size)) < array_size*array_size) {
                            wrap_around_z = false;
                            current = current + (array_size -1)*array_size*array_size;
                            } else current -= array_size*array_size;
                            hops++;
                            delta_z++;
                            }
                    }
                    
                    if (delta_z != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }

                    // X 
                    if (delta_x > 0){
                            while ( delta_x != 0 ){ // +x
                            int t = current * (degree+1+2*Host_Num) + 1;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+1+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 1;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( ((current%(array_size*array_size*array_size)) % (array_size*array_size)) % array_size == array_size-1) { 
                            wrap_around_x = false;
                            current = current - (array_size -1);
                            } else current++; 
                            delta_x--;
                            hops++;
                            }
                    } else if (delta_x < 0){
                            while ( delta_x != 0 ){ // -x
                            int t = current * (degree+1+2*Host_Num) + 2;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+2+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 2;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( ((current%(array_size*array_size*array_size)) % (array_size*array_size)) % array_size == 0 ) { 
                            wrap_around_x = false;
                            current = current + (array_size - 1);
                            } else current--;
                            hops++;
                            delta_x++;
                            }
                    }
                    
                    if (delta_x != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }

                    // Y 
                    if (delta_y > 0){
                            while ( delta_y != 0 ){ // -y
                            int t = current * (degree+1+2*Host_Num) + 3;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+3+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 3;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( ((current%(array_size*array_size*array_size)) % (array_size*array_size)) >= array_size*(array_size-1) ){ 
                            wrap_around_y = false;
                            current = current - array_size*(array_size -1);
                            } else current += array_size;
                            hops++;
                            delta_y--;
                            }
                    } else if (delta_y < 0){
                            while ( delta_y != 0 ){ // +y
                            int t = current * (degree+1+2*Host_Num) + 4;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+4+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 4;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( ((current%(array_size*array_size*array_size)) % (array_size*array_size)) < array_size ) {
                            wrap_around_y = false;
                            current = current + array_size*(array_size -1);
                            } else current -= array_size;
                            hops++;
                            delta_y++;
                            }
                    }
                    
                    if ( delta_x != 0 || delta_y != 0 || delta_z != 0 || delta_a != 0){ 
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }         
            }
            
            if (dimension == 3){ //3D routing
                    if (delta_z > 0){
                            while ( delta_z != 0 ){ // -z
                            int t = current * (degree+1+2*Host_Num) + 5;
                            // int t = (wrap_around_z) ? Vch*current*(degree+1+2*Host_Num)+5+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 5;
                            Crossing_Paths[t].pair_index.push_back(ct);  // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == array_size-1) { 
                            if ( current >= array_size*array_size*(array_size-1)) {
                            wrap_around_z = false;
                            current = current - (array_size -1)*array_size*array_size;
                            } else current += array_size*array_size; 
                            delta_z--;
                            hops++;
                            }
                    } else if (delta_z < 0){
                            while ( delta_z != 0 ){ // +z
                            int t = current * (degree+1+2*Host_Num) + 6;
                            // int t = (wrap_around_z) ? Vch*current*(degree+1+2*Host_Num)+6+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 6;
                            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            //if ( current % (array_size*array_size) == 0 ) { 
                            if ( current < array_size*array_size) {
                            wrap_around_z = false;
                            current = current + (array_size -1)*array_size*array_size;
                            } else current -= array_size*array_size;
                            hops++;
                            delta_z++;
                            }
                    }
                    
                    if (delta_z != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }

                    // X
                    if (delta_x > 0){
                            while ( delta_x != 0 ){ // +x
                            int t = current * (degree+1+2*Host_Num) + 1;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+1+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 1;
                            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( (current % (array_size*array_size)) % array_size == array_size-1) { 
                            wrap_around_x = false;
                            current = current - (array_size -1);
                            } else current++; 
                            delta_x--;
                            hops++;
                            }
                    } else if (delta_x < 0){
                            while ( delta_x != 0 ){ // -x
                            int t = current * (degree+1+2*Host_Num) + 2;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+2+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 2;
                            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( (current % (array_size*array_size)) % array_size == 0 ) { 
                            wrap_around_x = false;
                            current = current + (array_size - 1);
                            } else current--;
                            hops++;
                            delta_x++;
                            }
                    }
                    
                    // check X routing is finished 
                    if (delta_x != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }

                    // Y 
                    if (delta_y > 0){
                            while ( delta_y != 0 ){ // -y
                            int t = current * (degree+1+2*Host_Num) + 3;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+3+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 3;
                            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( (current % (array_size*array_size)) >= array_size*(array_size-1) ){ 
                            wrap_around_y = false;
                            current = current - array_size*(array_size -1);
                            } else current += array_size;
                            hops++;
                            delta_y--;
                            }
                    } else if (delta_y < 0){
                            while ( delta_y != 0 ){ // +y
                            int t = current * (degree+1+2*Host_Num) + 4;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+4+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 4;
                            Crossing_Paths[t].pair_index.push_back(ct);  // channel <-- node pair ID
                            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( (current % (array_size*array_size)) < array_size ) { 
                            wrap_around_y = false;
                            current = current + array_size*(array_size -1);
                            } else current -= array_size;
                            hops++;
                            delta_y++;
                            }
                    }
                    
                    // check if X,Y,Z routing are finished 
                    if ( delta_x != 0 || delta_y != 0 || delta_z != 0){ 
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }        
            }

            if (dimension == 2 || dimension == 1){ //2D routing, 1D routing (delta_y=0)
                    // X 
                    if (delta_x > 0){
                            while ( delta_x != 0 ){ // +x
                            int t = current * (degree+1+2*Host_Num) + 1;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+1+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 1;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( current % array_size == array_size-1) {
                            wrap_around_x = false;
                            current = current - (array_size -1);
                            } else current++; 
                            delta_x--;
                            hops++;
                            }
                    } else if (delta_x < 0){
                            while ( delta_x != 0 ){ // -x
                            int t = current * (degree+1+2*Host_Num) + 2;
                            // int t = (wrap_around_x) ? Vch*current*(degree+1+2*Host_Num)+2+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 2;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( current % array_size == 0 ) {
                            wrap_around_x = false;
                            current = current + (array_size - 1 );
                            } else current--;
                            hops++;
                            delta_x++;
                            }
                    }
                    
                    if (delta_x != 0){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }

                    // Y
                    if (delta_y > 0){
                            while ( delta_y != 0 ){ // -y
                            int t = current * (degree+1+2*Host_Num) + 3;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+3+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 3;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( current >= array_size*(array_size-1) ){
                            wrap_around_y = false;
                            current = current - array_size*(array_size -1);
                            } else current += array_size;
                            hops++;
                            delta_y--;
                            }
                    } else if (delta_y < 0){
                            while ( delta_y != 0 ){ // +y
                            int t = current * (degree+1+2*Host_Num) + 4;
                            // int t = (wrap_around_y) ? Vch*current*(degree+1+2*Host_Num)+4+(degree+1+2*Host_Num) :
                            // Vch * current * (degree+1+2*Host_Num) + 4;
                            Crossing_Paths[t].pair_index.push_back(ct); 
                            job.src_dst_pairs_m[i].channels.push_back(t);

                            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                            if ( current < array_size ) {
                            wrap_around_y = false;
                            current = current + array_size*(array_size -1);
                            } else current -= array_size;
                            hops++;
                            delta_y++;
                            }
                    }
                    
                    if ( delta_x != 0 || delta_y != 0 ){
                            cerr << "Routing Error " << endl;
                            exit (1);
                    }        
            }      

            // switch->host 
            t = dst*(degree+1+2*Host_Num)+degree+1+Host_Num+h_dst%Host_Num;
            // t = (src%2==1 && Topology==1) ? 
            //                 Vch*dst*(degree+1+2*Host_Num)+(degree+1+2*Host_Num)+degree+1+Host_Num+h_dst%Host_Num
            //                 : Vch*dst*(degree+1+2*Host_Num)+degree+1+Host_Num+h_dst%Host_Num;
            Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
            job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID   

            Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
            job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

            //ct++;	

        }
   }
   
   if (Topology == 2){ // fat-tree

    for (int i=0; i<job.src_dst_pairs_m.size(); i++){
        h_src = job.src_dst_pairs_m[i].h_src;
        h_dst = job.src_dst_pairs_m[i].h_dst;
        src = job.src_dst_pairs_m[i].src;
        dst = job.src_dst_pairs_m[i].dst;
        ct = job.src_dst_pairs_m[i].pair_id;

        int current = h_src;

        //#######################//
        // switch port:0 UP, 1 DOWN1(or localhost), 2 DOWN2, 3 DOWN3, 4 DOWN4 //
        //#######################//
        
        // channel --> switch ID + output port
        //Pair tmp_pair(src,dst,h_src,h_dst);  
        //pairs.push_back(tmp_pair);
        Crossing_Paths[current*PORT].pair_index.push_back(ct);
        job.src_dst_pairs_m[i].channels.push_back(current*PORT);        
        //hops++;

        //pairs[ct].pair_id = ct; 
        
        current = node_num + current/Host_Num;

        while ( current != dst ){ 
                int t;
                // root switch
                if ( current == node_num + node_num/Host_Num + node_num/pow(Host_Num,2)){
                t = current * PORT + dst/(int)pow(Host_Num,2)+1;
                current = current - Host_Num + dst/(int)pow(Host_Num,2);
                // middle layer switch
                } else if ( current >= node_num + node_num/Host_Num){
                if ( current-node_num-node_num/Host_Num != dst/(int)pow(Host_Num,2)){
                t = current * PORT + 0;
                current = node_num + node_num/Host_Num + node_num/(int)pow(Host_Num,2);
                } else {
                t = current * PORT + (dst/Host_Num)%Host_Num + 1;
                current = node_num + dst/Host_Num;
                }
                // low layer switch
                } else if ( current >= node_num){
                if ( current-node_num != dst/Host_Num){
                t = current * PORT + 0;
                current = node_num + current/Host_Num;
                } else {
                t = current * PORT + dst%Host_Num+1;
                current = dst;
                }
                }
                else { //host->switch
	        }      
                Crossing_Paths[t].pair_index.push_back(ct);
                job.src_dst_pairs_m[i].channels.push_back(t);

                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                hops++;
        }

        job.src_dst_pairs_m[i].hops = hops-before_hops;  
        before_hops = hops;

        if ( current != dst ){
                cerr << "Routing Error " << endl;
                exit (1);
        }
        //ct++;

    }	
   }
   
   
   if (Topology == 3) //fully-connected
   {
    for (int i=0; i<job.src_dst_pairs_m.size(); i++){
        h_src = job.src_dst_pairs_m[i].h_src;
        h_dst = job.src_dst_pairs_m[i].h_dst;
        src = job.src_dst_pairs_m[i].src;
        dst = job.src_dst_pairs_m[i].dst;
        ct = job.src_dst_pairs_m[i].pair_id;

        //#######################//
        // switch port <-- destination switch ID
        // localhost port = switch ID
        //#######################//

        // channel <-- node pair ID, node pair <-- channel ID 
        //Pair tmp_pair(src,dst,h_src,h_dst);  
        //pairs.push_back(tmp_pair);

        // localhost(h_src) --> src
        //int t = src*(degree+1+2*Host_Num)+degree+1+h_src%Host_Num;
        //Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        //pairs[ct].channels.push_back(t);  // node pair <-- channel ID     
        //pairs[ct].pair_id = ct; 
        job.src_dst_pairs_m[i].hops = 1;

        // src --> dst
        int t = src * (degree+1+2*Host_Num) + dst; // output port = destination switch ID
        Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        job.src_dst_pairs_m[i].channels.push_back(t);  // node pair <-- channel ID 

        Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
        job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

        hops++;      

        // dst --> localhost(h_dst)
        //t = dst*(degree+1+2*Host_Num)+degree+1+Host_Num+h_dst%Host_Num;
        //Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        //Crossing_Paths[dst*(degree+1+2*Host_Num)+dst].pair_index.push_back(ct); 
        //pairs[ct].channels.push_back(t); // node pair <-- channel ID  
        //t = dst*(degree+1+2*Host_Num)+dst;   
        //pairs[ct].channels.push_back(t); // output port <-- destination switch ID

        //ct++;
    }	
   }

      
   if (Topology == 4){ // full mesh connected circles (FCC)
        
        // switch topo initiation
        // switch n <--> switch 3-n, port p <--> port 11-p
        // 0:not used, 1:left, 2:right, 3-8:inter-group
        for(int g=0; g<groups; g++){
                for(int s=0; s<group_switch_num; s++){
                        int sw = g*group_switch_num+s;
                        if(s==0) Switch_Topo[sw*(degree+1+2*Host_Num)+1] = sw+(group_switch_num-1);
                        else Switch_Topo[sw*(degree+1+2*Host_Num)+1] = sw-1;
                        if(s==group_switch_num-1) Switch_Topo[sw*(degree+1+2*Host_Num)+2] = sw-(group_switch_num-1);
                        else Switch_Topo[sw*(degree+1+2*Host_Num)+2] = sw+1;
                        for(int j=3; j<=degree; j++){
                                Switch_Topo[sw*(degree+1+2*Host_Num)+j] = ((g+s*inter_group+j-3+1)%groups)*group_switch_num+(group_switch_num-1-s);
                        }
                }
        }

    for (int i=0; i<job.src_dst_pairs_m.size(); i++){
        h_src = job.src_dst_pairs_m[i].h_src;
        h_dst = job.src_dst_pairs_m[i].h_dst;
        src = job.src_dst_pairs_m[i].src;
        dst = job.src_dst_pairs_m[i].dst;
        ct = job.src_dst_pairs_m[i].pair_id;

        int current = src;   

        //#######################//
        // switch port: 0 localhost, 1 left, 2 right, 3-8 inter-group
        // localhost port = 0
        //#######################//

        // channel <-- node pair ID, node pair <-- channel ID 
        //Pair tmp_pair(src,dst,h_src,h_dst);  
        //pairs.push_back(tmp_pair);

        // localhost(h_src) --> src
        int t = src*(degree+1+2*Host_Num)+degree+1+h_src%Host_Num;
        Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        job.src_dst_pairs_m[i].channels.push_back(t);  // node pair <-- channel ID  

        Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
        job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);        

        //pairs[ct].pair_id = ct; 

        // src --> dst   
        // group #
        int src_group = src/group_switch_num;
        int dst_group = dst/group_switch_num;
        // switch offset in a group
        int src_offset = src%group_switch_num;
        int dst_offset = dst%group_switch_num;

        int diff_group = dst_group-src_group;
        // gateway for inter-group
        int src_gw_offset = -1;
        int dst_gw_offset = -1;
        int src_gw = -1;
        int dst_gw = -1;
        int src_gw_port = -1;
        int dst_gw_port = -1;
        if (diff_group==0){ //intra-group routing
                if(dst_offset-src_offset==0) continue;
                else if(dst_offset-src_offset>group_switch_num/2 || (src_offset-dst_offset>0 && src_offset-dst_offset<group_switch_num/2)){
                        while ( current != dst ){
                                t = current*(degree+1+2*Host_Num)+1; //backward
                                current = current - 1;
                                if(current<src_group*group_switch_num) current = current+group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }
                else{
                        while ( current != dst ){
                                t = current*(degree+1+2*Host_Num)+2; //forward
                                current = current + 1;
                                if(current>=(src_group+1)*group_switch_num) current = current-group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }
        }
        else { //intra-src-group routing + inter-group routing + intra-dst-group routing
                if(diff_group<0) diff_group+=groups; 
                src_gw_offset = (diff_group-1)/inter_group;
                dst_gw_offset = (group_switch_num-1)-src_gw_offset;
                src_gw = src_group*group_switch_num+src_gw_offset;
                dst_gw = dst_group*group_switch_num+dst_gw_offset;
                src_gw_port = (diff_group-1)%inter_group+3;
                dst_gw_port = (3+degree)-src_gw_port;

                //intra-src-group routing
                if(src_gw_offset-src_offset==0) ;
                else if(src_gw_offset-src_offset>group_switch_num/2 || (src_offset-src_gw_offset>0 && src_offset-src_gw_offset<group_switch_num/2)){
                        while ( current != src_gw ){
                                t = current*(degree+1+2*Host_Num)+1; //backward
                                current = current - 1;
                                if(current<src_group*group_switch_num) current = current+group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }
                else{
                        while ( current != src_gw ){
                                t = current*(degree+1+2*Host_Num)+2; //forward
                                current = current + 1;
                                if(current>=(src_group+1)*group_switch_num) current = current-group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }

                //inter-src-dst-group routing
                t = current*(degree+1+2*Host_Num)+src_gw_port;
                current = dst_gw;
                Crossing_Paths[t].pair_index.push_back(ct);
                job.src_dst_pairs_m[i].channels.push_back(t);

                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                hops++;

                //intra-dst-group routing
                if(dst_offset-dst_gw_offset==0) ;
                else if(dst_offset-dst_gw_offset>group_switch_num/2 || (dst_gw_offset-dst_offset>0 && dst_gw_offset-dst_offset<group_switch_num/2)){
                        while ( current != dst ){
                                t = current*(degree+1+2*Host_Num)+1; //backward
                                current = current - 1;
                                if(current<dst_group*group_switch_num) current = current+group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }
                else{
                        while ( current != dst ){
                                t = current*(degree+1+2*Host_Num)+2; //forward
                                current = current + 1;
                                if(current>=(dst_group+1)*group_switch_num) current = current-group_switch_num;
                                Crossing_Paths[t].pair_index.push_back(ct);
                                job.src_dst_pairs_m[i].channels.push_back(t);

                                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);

                                hops++;
                        }
                }                
        }

        // dst --> localhost(h_dst)
        t = dst*(degree+1+2*Host_Num)+degree+1+Host_Num+h_dst%Host_Num;
        Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID  

        Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
        job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);        

        job.src_dst_pairs_m[i].hops = hops-before_hops;  
        before_hops = hops;

        if ( current != dst ){
                cerr << "Routing Error " << endl;
                exit (1);
        }
        //ct++;
    }
   }

   if (Topology == 5){ // topology file
        
        // switch topo initiation (sw-port), store ports described in topology files
        // 0 - (switch_num-1) --> port (except itself); switch_num, (switch_num+1) --> not used
        for (int i=0; i<Switch_Topo.size(); i++){
                Switch_Topo[i] = -1;
        }
        for (int i=0; i<switch_num; i++){
                for (int j=0; j<topo_file.size(); j=j+2){
                        if (topo_sws_uni[i] == topo_file[j]){
                                int connect_sw = -1;
                                int connect_port = -1;
                                if (j%4 == 0){
                                        connect_sw = topo_file[j+2];
                                        connect_port = topo_file[j+1];
                                }
                                else if (j%4 == 2){
                                        connect_sw = topo_file[j-2];
                                        connect_port = topo_file[j+1];
                                }
                                for (int k=0; k<switch_num; k++){
                                        if (topo_sws_uni[k] == connect_sw){
                                                Switch_Topo[i*((switch_num-1)+1+2*Host_Num)+k] = connect_port;
                                        }
                                }
                        }
                }
        }

        // Number of vertices in the graph 
        int V = switch_num;

        // graph initialization
        vector<int> graph(V*V);
        for (int i=0; i<graph.size(); i++){
                if (Switch_Topo[(i/V)*((switch_num-1)+1+2*Host_Num)+i%V] == -1){
                        graph[i] = 0;
                }
                else {
                        graph[i] = 1;
                }
        }

        // dijkstra for each pair, store intermediate and dst sws (not including src sw)
        vector< vector<int> > pair_path(V*V);
        for (int i=0; i<V; i++){
                 dijkstra(V, graph, i, pair_path); 
        } 

    for (int i=0; i<job.src_dst_pairs_m.size(); i++){
        h_src = job.src_dst_pairs_m[i].h_src;
        h_dst = job.src_dst_pairs_m[i].h_dst;
        src = job.src_dst_pairs_m[i].src;
        dst = job.src_dst_pairs_m[i].dst;
        ct = job.src_dst_pairs_m[i].pair_id;

        // channel <-- node pair ID, node pair <-- channel ID 
        //Pair tmp_pair(src,dst,h_src,h_dst);  
        //pairs.push_back(tmp_pair);

        //pair path
        int src_index = -1;
        int dst_index = -1;
        for (int i = 0; i < switch_num; i++){
                if (topo_sws_uni[i] == src){
                        src_index = i;
                }
                if (topo_sws_uni[i] == dst){
                        dst_index = i;
                } 
                if (src_index != -1 && dst_index != -1) break;               
        }
        if (src_index == -1 || dst_index == -1) cout << "Error: src/dst number is wrong" << endl;

        // localhost(h_src) --> src
        int t = src_index*((switch_num-1)+1+2*Host_Num)+(switch_num-1)+1+h_src%Host_Num;
        Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        job.src_dst_pairs_m[i].channels.push_back(t);  // node pair <-- channel ID   

        Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
        job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);        

        //pairs[ct].pair_id = ct; 
        job.src_dst_pairs_m[i].hops = pair_path[src_index*V+dst_index].size()+1;
        hops += job.src_dst_pairs_m[i].hops;

        // src --> dst
        for (int i=0; i<pair_path[src_index*V+dst_index].size(); i++){
                if (i==0){
                        t = src_index*((switch_num-1)+1+2*Host_Num)+pair_path[src_index*V+dst_index][i];
                }
                else {
                        t = pair_path[src_index*V+dst_index][i-1]*((switch_num-1)+1+2*Host_Num)+pair_path[src_index*V+dst_index][i];
                }
                Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
                job.src_dst_pairs_m[i].channels.push_back(t);  // node pair <-- channel ID 

                Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
                job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);                

        }

        // dst --> localhost(h_dst)
        t = dst_index*((switch_num-1)+1+2*Host_Num)+(switch_num-1)+1+Host_Num+h_dst%Host_Num;
        Crossing_Paths[t].pair_index.push_back(ct); // channel <-- node pair ID
        job.src_dst_pairs_m[i].channels.push_back(t); // node pair <-- channel ID  

        Crossing_Paths[t].flow_index.push_back(job.job_id*1000+job.src_dst_pairs_m[i].flow_id);
        job.flows[job.src_dst_pairs_m[i].flow_id].channels.push_back(t);        

        //ct++;
    } 
   }
}

void calc_slots(int ports, vector<Cross_Paths> & Crossing_Paths, Job & job, vector<Job> & all_jobs){

    // delete duplicate flow id in channel
    for (int i = 0; i < ports; i++){
            vector<Cross_Paths>::iterator elem = Crossing_Paths.begin()+i;
            sort(elem->flow_index.begin(), elem->flow_index.end());
            elem->flow_index.erase(unique(elem->flow_index.begin(), elem->flow_index.end()), elem->flow_index.end());
    }

    // delete duplicate channel id in flow
    for (int i = 0; i < job.flows.size(); i++){
            vector<Flow>::iterator elem = job.flows.begin()+i;
            sort(elem->channels.begin(), elem->channels.end());
            elem->channels.erase(unique(elem->channels.begin(), elem->channels.end()), elem->channels.end());           
    }

   int max_cp = 0;
   int max_cp_dst = 0; // number of dst-based renewable labels
   // calculate number of dst-based renewable label
   int max_cp_dst_t = 0;

    for (int j = 0; j < ports; j++ ){ 
            vector<Cross_Paths>::iterator elem = Crossing_Paths.begin()+j;
            if (elem->flow_index.size() > max_cp) max_cp = elem->flow_index.size();
            unsigned int p_ct = 0;
            while ( p_ct < elem->flow_index.size() ){
                    int u = elem->flow_index[p_ct]; 
                    bool is_duplicate = false; //check if there is a same destination
                    unsigned int p_ct_t = 0;
                    while ( p_ct_t < p_ct ){
                            int v = elem->flow_index[p_ct_t]; 
                            for (int m = 0; m < all_jobs[u/1000].src_dst_pairs_m.size(); m++){
                                if (all_jobs[u/1000].src_dst_pairs_m[m].flow_id == u%1000){
                                    for (int n = 0; n < all_jobs[v/1000].src_dst_pairs_m.size(); n++){
                                        if (all_jobs[v/1000].src_dst_pairs_m[n].flow_id == v%1000){
                                            if (all_jobs[u/1000].src_dst_pairs_m[m].h_dst == all_jobs[v/1000].src_dst_pairs_m[n].h_dst){
                                                    is_duplicate = true;
                                                    break;
                                            }
                                        }
                                    }
                                    if (is_duplicate == true) break;
                                }
                            }
                            if (is_duplicate == true) break;
                            p_ct_t++;
                    }
                    if (!is_duplicate) max_cp_dst_t++;
                    p_ct++;
            }
            if (max_cp_dst_t > max_cp_dst) max_cp_dst = max_cp_dst_t;
            max_cp_dst_t = 0;
    }
    cout << " === Max. number of slots (w/o update) [FLOW] ===" << endl << max_cp_dst << endl;	
    cout << " === Max. number of slots (w/ update) [FLOW] ===" << endl << max_cp << endl;
}

void alloc_slot(Job & job, vector<Cross_Paths> & Crossing_Paths, bool path_based, int & max_id){
    for (int i=0; i<job.flows.size(); i++){
        
        // local IDs are assigned
        unsigned int path_ct = 0; 
        while ( path_ct < job.flows.size() ){   
            // check if IDs are assigned
            bool valid = true;

            for (int n=0; n<job.src_dst_pairs_m.size(); n++){
                if (job.src_dst_pairs_m[n].flow_id == job.flows[path_ct].id){
                    if (job.src_dst_pairs_m[n].Valid == false){
                        valid = false;
                        break;
                    }
                }
            }

            if ( valid == true ) {path_ct++; continue;}
            // ID is assigned from 0
            int id_tmp = 0;
            bool NG_ID = false;
                    
            NEXT_ID_FLOW:
                    // ID is used or not
                    unsigned int s_ct = 0; // channel
                    while ( s_ct < job.flows[path_ct].channels.size() && !NG_ID ){
                        int i = job.flows[path_ct].channels[s_ct];
                        vector<int>::iterator find_ptr;
                        find_ptr = find ( Crossing_Paths[i].assigned_list.begin(), Crossing_Paths[i].assigned_list.end(), id_tmp);
                        if ( path_based && find_ptr != Crossing_Paths[i].assigned_list.end()) NG_ID = true;
                        if (!path_based && find_ptr != Crossing_Paths[i].assigned_list.end()) {
                            int tmp = 0;
                            while (*find_ptr != Crossing_Paths[i].assigned_list[tmp]) {tmp++;}
                            NG_ID = true; 
                            for (int j=0; j<job.src_dst_pairs_m.size(); j++){
                                if (job.src_dst_pairs_m[j].flow_id == job.flows[path_ct].id){
                                    if (job.src_dst_pairs_m[j].h_dst == Crossing_Paths[i].assigned_dst_list[tmp]){
                                            NG_ID = false;
                                            break;
                                    }
                                }
                            }
                        }
                        s_ct++;
                    }
                    if (NG_ID){
                    id_tmp++; NG_ID = false; goto NEXT_ID_FLOW;
                    }
                    job.flows[path_ct].ID = id_tmp;

                    unsigned int a_ct = 0;
                    while ( a_ct < job.flows[path_ct].channels.size() ){
                            int j = job.flows[path_ct].channels[a_ct];
                            Crossing_Paths[j].assigned_list.push_back(id_tmp); 
                            for (int n=0; n<job.src_dst_pairs_m.size(); n++){
                                if (job.src_dst_pairs_m[n].flow_id == job.flows[path_ct].id){
                                    Crossing_Paths[j].assigned_dst_list.push_back(job.src_dst_pairs_m[n].h_dst);
                                }
                            }                                                               	
                            a_ct++; 
                    }

                    for (int n=0; n<job.src_dst_pairs_m.size(); n++){
                        if (job.src_dst_pairs_m[n].flow_id == job.flows[path_ct].id){
                            job.src_dst_pairs_m[n].Valid = true;
                        }
                    }
                                  	                       	    
                    if (max_id <= id_tmp) max_id = id_tmp + 1; 

                    path_ct++;
        }
        //elem->Valid = true;

    }
}

void dispatch_random(vector<int> & sys, vector<Job> & queue, vector<Job> & all_jobs, int topo, int switch_num, int hosts, int degree, vector<Cross_Paths> & Crossing_Paths, int dimension, int array_size, int & hops, int & before_hops, int PORT, int node_num,
                        int groups, int group_switch_num, vector<int> & Switch_Topo, int inter_group, vector<int> topo_file, vector<int> topo_sws_uni,
                        int & src, int & dst, int & h_src, int & h_dst, int & ct, int ports, bool path_based, int & max_id){
        vector<int> ava_nodes;
        for (int i=0; i<sys.size(); i++){
                if(sys[i] == 0){ //available
                       ava_nodes.push_back(i); 
                       if (queue[0].num_nodes == ava_nodes.size()){
                               //cout << "Job " << queue[0].job_id << " is dispatched" << endl;
                               printf("Job %d is dispatched\n", queue[0].job_id); 
                               //for (int j=0; j<ava_nodes.size(); j++) cout << ava_nodes[i] << " " << endl;
                               for (int j=0; j<all_jobs.size(); j++){
                                       if(all_jobs[j].job_id == queue[0].job_id){
                                                all_jobs[j].time_dispatch_r = time(NULL);
                                                for (int n=0; n<ava_nodes.size(); n++){
                                                        sys[ava_nodes[n]] = 1; //occupied
                                                        all_jobs[j].nodes.push_back(ava_nodes[n]);
                                                } 
                                                all_jobs[j].src_dst_pairs_m = src_dst_pairs_map(all_jobs[j].src_dst_pairs, all_jobs[j].nodes, topo, hosts);
                                                
                                                //todo
                                                route_job(topo, all_jobs[j], degree, switch_num, hosts, Crossing_Paths, dimension, array_size, hops, before_hops, PORT, node_num, 
                                                            groups, group_switch_num, Switch_Topo, inter_group, topo_file, topo_sws_uni, 
                                                            src, dst, h_src, h_dst, ct);
                                                calc_slots(ports, Crossing_Paths, all_jobs[j], all_jobs);
                                                alloc_slot(all_jobs[j], Crossing_Paths, path_based, max_id);

                                                thread rn{release_nodes, all_jobs[j], ref(sys), ref(Crossing_Paths)};
                                                rn.detach();  
                                                queue.erase(queue.begin()); 
                                                return;          
                                       }
                               }       
                       }
                }
        }
}

int main(int argc, char *argv[])
{
   static int array_size = 4;
   static int Host_Num = 1;
   static int Topology = 0; 

   static int queue_policy = 1;

        // int c;
        // static char* topology_file; // topology file

        // while((c = getopt(argc, argv, "t:")) != -1) {
        // switch (c) { 
        //         case 't': //topology file (Topology = 5)
        //                 topology_file = optarg;
        //                 break;               
        //         default:
        //                 //usage(argv[0]);
        //                 cout << " This option is not supported. " << endl;
        //                 return EXIT_FAILURE;
        //         }
        // }

        // vector<int> topo_file;
        // vector<int> topo_sws_dup; // duplicate
        // vector<int> topo_sws_uni; // unique
        // ifstream infile_feat(topology_file);
        // string line_data;
        // int data;
        // while (!infile_feat.eof()){
        //         getline(infile_feat, line_data);
        //         stringstream stringin(line_data);
        //         int column = 0;
        //         while (stringin >> data){
        //                 if (column == 0){ // src_sw
        //                         topo_file.push_back(data);
        //                         topo_sws_dup.push_back(data);
        //                         topo_sws_uni.push_back(data);
        //                         column++;
        //                 }
        //                 else if (column == 1){ // src_port
        //                         topo_file.push_back(data);    
        //                         column++;                           
        //                 }
        //                 else if (column == 2){ // dst_sw
        //                         topo_file.push_back(data); 
        //                         topo_sws_dup.push_back(data);
        //                         topo_sws_uni.push_back(data);  
        //                         column++;                                                           
        //                 }
        //                 else if (column == 3){ // dst_port
        //                         topo_file.push_back(data);   
        //                         column++;                            
        //                 }                                                
        //         }
        // }
        // infile_feat.close();

        // sort(topo_sws_uni.begin(), topo_sws_uni.end());
        // topo_sws_uni.erase(unique(topo_sws_uni.begin(), topo_sws_uni.end()), topo_sws_uni.end());
        // int switch_num = topo_sws_uni.size();
        // int degree = 0;
        // int cnt;
        // for (int i=0; i<switch_num; i++){
        //         cnt = count(topo_sws_dup.begin(), topo_sws_dup.end(), topo_sws_uni[i]);
        //         if (cnt > degree) degree = cnt;
        // }

        // for (int i=0; i<topo_file.size(); i++) cout << topo_file[i] << " ";
        // cout << endl;
        // for (int i=0; i<topo_sws_uni.size(); i++) cout << topo_sws_uni[i] << " ";
        // cout << endl;   
        // for (int i=0; i<topo_sws_dup.size(); i++) cout << topo_sws_dup[i] << " ";    
        // cout << endl; 
        // cout << degree << endl;      
        // cout << switch_num << endl;    

        // static int Host_Num = 1;
        // int ports = ((switch_num-1)+1+2*Host_Num)*switch_num;  
        // vector<int> Switch_Topo(ports);

        // for (int i=0; i<Switch_Topo.size(); i++){
        //         Switch_Topo[i] = -1;
        // }
        // for (int i=0; i<switch_num; i++){
        //         for (int j=0; j<topo_file.size(); j=j+2){
        //                 if (topo_sws_uni[i] == topo_file[j]){
        //                         int connect_sw = -1;
        //                         int connect_port = -1;
        //                         if (j%4 == 0){
        //                                 connect_sw = topo_file[j+2];
        //                                 connect_port = topo_file[j+1];
        //                         }
        //                         else if (j%4 == 2){
        //                                 connect_sw = topo_file[j-2];
        //                                 connect_port = topo_file[j+1];
        //                         }
        //                         for (int k=0; k<switch_num; k++){
        //                                 if (topo_sws_uni[k] == connect_sw){
        //                                         Switch_Topo[i*((switch_num-1)+1+2*Host_Num)+k] = connect_port;
        //                                 }
        //                         }
        //                 }
        //         }
        // }

        // for (int i=0; i<Switch_Topo.size(); i++) {
        //         cout << Switch_Topo[i] << " ";    
        //         if (i%((switch_num-1)+1+2*Host_Num) == ((switch_num-1)+1+2*Host_Num)-1)
        //                 cout << endl;  
        // }

        // int V = switch_num;
        // vector<int> graph(V*V);
        // for (int i=0; i<graph.size(); i++){
        //         if (Switch_Topo[(i/V)*((switch_num-1)+1+2*Host_Num)+i%V] == -1){
        //                 graph[i] = 0;
        //         }
        //         else {
        //                 graph[i] = 1;
        //         }
        // }
        // for (int i = 0; i < graph.size(); i++){
        //         cout << graph[i] << " ";
        //         if (i%V == V-1) cout << endl;
        // }

        // vector< vector<int> > pair_path(V*V);
        // for (int i=0; i<V; i++){
        //          dijkstra(V, graph, i, pair_path); 
        // }  
        // //dijkstra(V, graph, 1, pair_path); 

        // for (int i=0; i<V*V; i++){
        //         cout << i << ": ";
        //         for (int j=0; j<pair_path[i].size(); j++){
        //                 cout << pair_path[i][j] << " ";
        //         }
        //         cout << endl;
        // } 

        // int a = -1, b = -1, c = -1;
        // while ( cin >> a){	
        //         cin >> b; 
        //         // getchar();
        //         // char x;
        //         // while((x=cin.get()) != '\n') 
        //         //     if(x != ' ' && x != '\t'){
        //         //         //c = x-48;
        //         //         cin >> c;
        //         //         break;
        //         //     }   
        //         char array[20];
        //         cin.getline(array, 20);
        //         stringstream stringin(array);
        //         int flowid;
        //         while (stringin >> flowid){c = flowid;}
        //         cout << a << " " << b << " " << c << endl;
        // }     

   static int degree = 4; 
   static int dimension = 2; 

   static int group_switch_num = 4; 
   static int inter_group = 6; 
   static int intra_group = 2; 

   int src = -1, dst = -1, h_src = -1, h_dst = -1;

   int switch_num = 10;
   static int node_num = 10;
   static int PORT = Host_Num + 1; //fat-tree

   static int groups = inter_group*group_switch_num+1; //fcc 25

   vector<int> sys(switch_num);

   int temp_time_submit, temp_time_run;
   int temp_num_nodes, temp_pair_src, temp_pair_dst, temp_flow_id, temp_job_id;
   int pre_job_id = -1;
   int pre_flow_id = -1;

   vector<Job> all_jobs;
   vector<Job> queue;

   vector<int> topo_file;
   vector<int> topo_sws_dup; // duplicate
   vector<int> topo_sws_uni; // unique

   static int ports = 7*switch_num;
   vector<Cross_Paths> Crossing_Paths(ports);

   vector<int> Switch_Topo(ports);

   static bool path_based = false; 

   int max_id = 0;

   int ct = 0;
   int hops = 0; 
   int before_hops = 0;    

   int pairID = 0;
   while (cin >> temp_time_submit >> temp_time_run >> temp_num_nodes >> temp_pair_src >> temp_pair_dst >> temp_flow_id >> temp_job_id){
        if (temp_job_id != pre_job_id){
                Job j;
                j.time_submit = temp_time_submit;
                j.time_run = temp_time_run;
                j.num_nodes = temp_num_nodes;

                // j.src_dst_pair.push_back(temp_pair_src);
                // j.src_dst_pair.push_back(temp_pair_dst);
                // j.src_dst_pair.push_back(temp_flow_id);
                // j.src_dst_pair.push_back(pairID);
                // pairID++;
                // j.src_dst_pair.push_back(-1); // slot #
                int h_src = temp_pair_src;
                int h_dst = temp_pair_dst;
                int src = h_src/Host_Num;
                int dst = h_dst/Host_Num;
                if(Topology == 2){ // fat-tree
                     src = h_src;
                     dst = h_dst;   
                }
                Pair pair(src,dst,h_src,h_dst);
                pair.pair_id = pairID;
                pair.flow_id = temp_flow_id;
                j.src_dst_pairs.push_back(pair);

                Flow f;
                f.id = temp_flow_id;
                f.ID = -1;
                f.pairs_id.push_back(pairID);
                j.flows.push_back(f);

                j.job_id = temp_job_id;
                j.time_submit_r = -1;
                j.time_dispatch_r = -1;
                all_jobs.push_back(j);

                pre_job_id = temp_job_id;
                pre_flow_id = temp_flow_id;

                pairID++;
        }
        else{
                //src0, dst0, flow id0, pairID0, slot #0, src1, dst1, flow id1, pairID1, slot #1, ...
                // all_jobs[all_jobs.size()-1].src_dst_pair.push_back(temp_pair_src);
                // all_jobs[all_jobs.size()-1].src_dst_pair.push_back(temp_pair_dst);
                // all_jobs[all_jobs.size()-1].src_dst_pair.push_back(temp_flow_id);
                // all_jobs[all_jobs.size()-1].src_dst_pair.push_back(pairID);
                // pairID++;
                // all_jobs[all_jobs.size()-1].src_dst_pair.push_back(-1); // slot #
                int h_src = temp_pair_src;
                int h_dst = temp_pair_dst;
                int src = h_src/Host_Num;
                int dst = h_dst/Host_Num;
                if(Topology == 2){ // fat-tree
                     src = h_src;
                     dst = h_dst;   
                }
                Pair pair(src,dst,h_src,h_dst);
                pair.pair_id = pairID;
                pair.flow_id = temp_flow_id;
                all_jobs[all_jobs.size()-1].src_dst_pairs.push_back(pair);  

                if (temp_flow_id != pre_flow_id){
                        Flow f;
                        f.id = temp_flow_id;
                        f.ID = -1;
                        f.pairs_id.push_back(pairID);
                        all_jobs[all_jobs.size()-1].flows.push_back(f);
                        pre_flow_id = temp_flow_id;                         
                }
                else {
                        all_jobs[all_jobs.size()-1].flows[all_jobs[all_jobs.size()-1].flows.size()-1].pairs_id.push_back(pairID);
                }

                pairID++;              
        }

   }

//    cout << all_jobs.size() << endl; 
//    for (int i=0; i<all_jobs.size(); i++){
//        cout << all_jobs[i].flows.size() << endl;
//    }  
//    for (int i=0; i<all_jobs[2].flows[2].pairs_id.size(); i++){
//        cout << all_jobs[2].flows[2].pairs_id[i] << endl;
//    }  

   //job submit thread
   bool all_submitted = false;
   thread sj{submit_jobs, ref(all_jobs), ref(queue), queue_policy, ref(all_submitted)};
   sj.detach();   

   //check system state
   bool all_finished = false;
   thread cs{check_state, ref(sys), ref(queue), ref(all_submitted), ref(all_finished)};
   cs.detach();     

//    vector<int> b;
//    b.push_back(1);
//    b.push_back(3);
//    greeting(2, b);
//    thread t{greeting, 2, ref(b)};
//    t.join(); 
//    for (int i=0; i<b.size(); i++){
//        cout<< b[i] << endl;
//    }

    // while (true) {
    //     if (all_submitted){
    //         cout << all_jobs[2].time_submit << endl;
    //         cout << all_jobs[2].num_nodes << endl;
    //         cout << all_jobs[2].job_id << endl;
    //         cout << all_jobs[2].time_submit_r << endl;
    //         for (int i=0; i<queue.size(); i++){
    //             cout << queue[i].job_id << endl;
    //         }
    //         break;
    //     }
    // }

   while(true){
        if(queue.size() > 0){
                if(queue[0].num_nodes<1 || queue[0].num_nodes>switch_num || queue[0].time_run<0){
                        cout << "Job " << queue[0].job_id << "can not be dispatched due to errorous request" << endl;
                        queue.erase(queue.begin());
                        //continue;
                }
                else{
                        dispatch_random(sys, queue, all_jobs, Topology, switch_num, Host_Num, degree, Crossing_Paths, dimension, array_size, hops, before_hops, PORT, node_num,
                        groups, group_switch_num, Switch_Topo, inter_group, topo_file, topo_sws_uni,
                        src, dst, h_src, h_dst, ct, ports, path_based, max_id);
                }
        }
        else if(all_submitted == true){
                cout << "All jobs have been dispatched" << endl;
                while (true) {
                    if (all_finished == true){
                        cout << "All jobs have been finished" << endl;
                        break;
                    }
                }
                break;  
        }
   }

    for (int i=0; i<all_jobs.size(); i++){
        for (int j=0; j<all_jobs[i].nodes.size(); j++)
        {
            cout << all_jobs[i].nodes[j] << "  ";
        }
        cout << endl;
    }   
    for (int i=0; i<all_jobs[2].src_dst_pairs.size(); i++){
        cout << all_jobs[2].src_dst_pairs[i].src << "  " << all_jobs[2].src_dst_pairs[i].dst << endl;
    }
    for (int i=0; i<all_jobs[2].src_dst_pairs_m.size(); i++){
        cout << all_jobs[2].src_dst_pairs_m[i].src << "  " << all_jobs[2].src_dst_pairs_m[i].dst << endl;
    }

   pthread_exit(NULL);    
   return 0;
}
