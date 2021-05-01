
#include "deform_graph.h"

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::BuildGraph(
    const TriMesh& mesh,
    pcl::PolygonMesh& in_graph, 
    int k_nn) 
{    
    MeshProcessing mp;
    CloudPtr graph_pts (new pcl::PointCloud<pcl::PointXYZRGBNormal>());

    graph_pts = mp.GetCloudFromPolygonMesh(in_graph);

    node_pos.clear(); 
    node_norm.clear();
    node_rots_mat.clear();
    node_rots.clear();
    node_trans.clear();
    node_neigh.clear();
    node_node_weights.clear();
    vert_node_weights.clear();

    node_pos.resize(graph_pts->size());
	node_norm.resize(graph_pts->size());
    
    for (size_t i=0; i<node_pos.size(); ++i) 
    {
        Eigen::Vector3d point, normal; 
        point.setZero(); normal.setZero();
        point(0) = graph_pts->points[i].x;
        point(1) = graph_pts->points[i].y;
        point(2) = graph_pts->points[i].z;
        normal(0) = graph_pts->points[i].normal_x;
        normal(1) = graph_pts->points[i].normal_y;
        normal(2) = graph_pts->points[i].normal_z;
        node_pos[i]  = point;
        node_norm[i] = normal;
    }
    
    // Init values
    node_rots_mat.resize(node_pos.size(),Eigen::Matrix3d::Identity());
    node_rots.resize(node_pos.size(), Eigen::Quaterniond::Identity());
    node_trans.resize(node_pos.size(), Eigen::Vector3d::Zero());

    // TODO: Add data structures for storing all node-node weights, dists etc. 
    // so that it can be done in a single (albeit long) operation
    
    this->SetMGraphStructure(in_graph);
    BuildNeighbours(); // knn build connectivity
    CalculateNodeNodeWeights();
}

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::BuildNeighbours() 
{
    node_neigh.resize(node_pos.size());
    MeshProcessing mp;

    Mesh he_mesh; 
    pcl::geometry::toHalfEdgeMesh(graph_structure_,he_mesh);

    // sanity
    if (node_pos.size() != he_mesh.sizeVertices())
    {
        std::string e = 
        "In BuildNeighbours: structure vertices size doesn't match graph nodes size";
        BOOST_LOG_TRIVIAL(error) << e;
        throw std::runtime_error(e);
    }
    
    // for each vertex we build a circulator and push back all 
    // connected neighbors
    for (size_t pt_idx=0; pt_idx<node_pos.size(); ++pt_idx) 
    {
        
        Mesh::VertexIndex vi(pt_idx);
        Mesh::VertexAroundVertexCirculator vavc = 
            he_mesh.getVertexAroundVertexCirculator(vi);

        const Mesh::VertexAroundVertexCirculator vavc_end = vavc;
        if (!vavc.isValid())
        {
            continue;
        }
        
        do 
        {
            node_neigh[pt_idx].push_back(vavc.getTargetIndex().get());
        } while (++vavc != vavc_end);

    }

    // TODO: revise this, technically the implementation calls that all 
    // nodes which influence a given point sharing a connecting edge. 
    // Current implementation will be severly lacking in connectivity 
}

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::CalculateNodeNodeWeights()
{
    if (this->node_neigh.empty())
    {
        BOOST_LOG_TRIVIAL(warning) 
            << "In CalculateNodeNodeWeights: Calculate weights called "
            << "before neighbour connectivity built. Building now...";
        BuildNeighbours();
    }

    // NodeNode weights won't vary much if the graph is uniformly distributed
    this->node_node_weights.clear();
    this->node_node_weights.resize(node_neigh.size());

    // node-node weights
    for (size_t node_idx = 0; node_idx < this->node_neigh.size(); node_idx++)
    {
        int num_neigh = node_neigh[node_idx].size(); // should be ~kDG_KNN
        std::vector<double> weights; 
        weights.resize(num_neigh);

        // get sqr distances to nn_nodes
        std::vector<double> dists_node; 
        dists_node.clear(); 
        dists_node.resize(num_neigh);
        Eigen::Vector3d pos = this->node_pos[node_idx];

        for(size_t neigh_idx = 0; neigh_idx < num_neigh; ++neigh_idx)
        {

            int nid = node_neigh[node_idx][neigh_idx];

            Eigen::Vector3d node_pos = this->node_pos[nid];
            double sq_dist = 0;
            sq_dist += (pos(0)-node_pos(0))*(pos(0)-node_pos(0));
            sq_dist += (pos(1)-node_pos(1))*(pos(1)-node_pos(1));
            sq_dist += (pos(2)-node_pos(2))*(pos(2)-node_pos(2));
            double dist = sqrt(sq_dist);

            dists_node.push_back(dist);
        }
        
        double weight_sum = 0;

        double d_max = 1.5 * this->node_sample_radius_;
        
        for(size_t neigh_idx = 0; neigh_idx < num_neigh; ++neigh_idx)
        {
            weights[neigh_idx] = utils::max(
                0.0, utils::cube(1.0 - utils::square(
                        dists_node[neigh_idx] / d_max)));         
            weight_sum += weights[neigh_idx];
        }        

        for(size_t neigh_idx = 0; neigh_idx < num_neigh; ++neigh_idx)
        {
            weights[neigh_idx] /= weight_sum;
        }

        // store
        this->node_node_weights[node_idx] = weights;
    }

}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> DeformGraph::CalculateVertNodeWeights(
    const std::vector<pcl::PointXYZRGBNormal>& verts,
    const std::vector<std::vector<int>>& vert_neigh)  
{

    std::vector<std::vector<double>> vert_node_weights;
    vert_node_weights.resize(verts.size());

    for (size_t con_idx = 0; con_idx < verts.size(); con_idx++)
    {
        std::vector<double> weights;
        weights.resize(vert_neigh[con_idx].size());

        // get sqr distances to nn_nodes
        std::vector<double> dists_node;
        Eigen::Vector3d pos; pos.setZero();
        pos(0) = verts[con_idx].x;
        pos(1) = verts[con_idx].y;
        pos(2) = verts[con_idx].z;

        // use euclidean dist
        for (size_t neigh_idx = 0; neigh_idx < weights.size(); neigh_idx++)
        {
            size_t nid = vert_neigh[con_idx][neigh_idx];
            Eigen::Vector3d node_pos = this->node_pos[nid];
            double sq_dist = 0;
            sq_dist += (pos(0)-node_pos(0))*(pos(0)-node_pos(0));
            sq_dist += (pos(1)-node_pos(1))*(pos(1)-node_pos(1));
            sq_dist += (pos(2)-node_pos(2))*(pos(2)-node_pos(2));
            double dist = sqrt(sq_dist);

            dists_node.push_back(dist);
        }
        double weight_sum = 0;

        // use a scalar to increase the influence radius as the sample radius
        // is proportional to the graph density
        double dmax_scaler = this->dmax_scaler_;

        // keep increasing sampling radius if we can't find a node
        while (weight_sum <= 0.5)
        {

            weight_sum = 0;
            for (size_t neigh_idx = 0; neigh_idx < weights.size(); neigh_idx++)
            {
                double d_max = dmax_scaler * this->node_sample_radius_;
                weights[neigh_idx] = utils::max(
                    0.0, utils::cube(1.0 - utils::square(
                            dists_node[neigh_idx] / d_max)));       
                weight_sum += weights[neigh_idx];
            }

            dmax_scaler += 0.05;
        }

        for (size_t neigh_idx = 0; neigh_idx < weights.size(); neigh_idx++)
        {
            weights[neigh_idx] /= weight_sum;
        }
              
        vert_node_weights[con_idx] = weights;

        // normalization check
        weight_sum = 0;
        for (size_t w = 0; w < weights.size(); w++)
        {
            weight_sum += weights[w];
        }
        // roughly equal to 1
        if (std::abs(weight_sum - 1) > 0.01)
        {
            std::string e 
                = "In DeformGraph::GetVertNodeWeights, weight normalization failed";
            BOOST_LOG_TRIVIAL(debug) << e <<", sum: "<<weight_sum;
            throw std::runtime_error(e);
        }
    }
    return vert_node_weights;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> DeformGraph::CalculateVertNodeWeights(
    const pcl::PolygonMesh& in_mesh,
    const std::vector<std::vector<int>>& vert_neigh)  
{
    MeshProcessing mp;
    pcl::PolygonMesh tmp_mesh (in_mesh);
    CloudPtr in_cloud = mp.GetCloudFromPolygonMesh(tmp_mesh);

    std::vector<pcl::PointXYZRGBNormal> verts;
    for (size_t pt = 0; pt < in_cloud->size(); pt++)
    {
        verts.push_back(in_cloud->points[pt]);
    }

    return CalculateVertNodeWeights(verts,vert_neigh);
}

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::UpdateGraph(
    const Eigen::VectorXd& X) 
{
	int num_nodes = node_pos.size();
	std::vector<Eigen::Matrix3d> rots(num_nodes); 
    std::vector<Eigen::Vector3d> trans(num_nodes);

	for (int node=0; node<num_nodes; ++node) 
    {
		for (int r=0; r<3; ++r) 
        {
			for (int c=0; c<3; ++c)
            {
				rots[node](r, c) = X[12*node+3*r+c];
            }
			trans[node][r] = X[12*node+9+r];
		}
		Eigen::Matrix3d r = rots[node].inverse(); 

        rots[node] = r;
		rots[node].transposeInPlace();

		// node_pos[node] = node_pos[node] + trans[node];
		node_norm[node] = rots[node]*node_norm[node];
	}
	for (int node=0; node<num_nodes; ++node) 
    {
		node_rots[node] = rots[node]*node_rots[node];
		node_trans[node] = rots[node]*node_trans[node] + trans[node];
	}
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>> DeformGraph::GetVertNodeNeighbours(
    pcl::PolygonMesh& in_mesh)
{

    std::vector<std::vector<int>> node_knn_lists;

    Matcher mt;
    MeshProcessing mp;
    CloudPtr mesh_cloud = mp.GetCloudFromPolygonMesh(in_mesh);
    CloudPtr query_cloud = mp.GetCloudFromPolygonMesh(this->graph_structure_);

    for (size_t vid = 0; vid < mesh_cloud->size(); vid++)
    {
        std::vector<int> nodes;

        try
        {
            nodes = GetVertNodeMatches(
                mesh_cloud->points[vid],
                this->graph_structure_,
                query_cloud);
        }
        catch(const std::exception& e)
        {
            BOOST_LOG_TRIVIAL(error) << "point: "<<vid;
            std::cin.get();

            // in case of fail copy neighbor nodes
            node_knn_lists.push_back(node_knn_lists.back());
            // throw std::runtime_error(e.what());
        }
        
        // TODO: Add flag to copy neighbor nodes
        if (nodes.size() == 0)            
        {
            std::string e = 
                "In DeformGraph::GetVertNodeNeighbours, no matches found";
            BOOST_LOG_TRIVIAL(warning) << e;
            nodes = node_knn_lists.back();
            // throw std::runtime_error(e);
        }

        node_knn_lists.push_back(nodes);
    }
    
    return node_knn_lists;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> DeformGraph::GetVertNodeMatches(
    const pcl::PointXYZRGBNormal& source, 
    const pcl::PolygonMesh& query,
    const CloudPtr& query_cloud)
{

    c_clock::time_point tic,toc;

    MeshProcessing mp;

    std::vector<int> matches;


    // sanity check 
    Eigen::Vector3d src_n{
        source.normal_x,
        source.normal_y,
        source.normal_z
    };

    if (src_n[0] == 0.0 && src_n[1] == 0.0 && src_n[2] == 0.0 )
    {
        BOOST_LOG_TRIVIAL(error) 
        << "In DeformGraph::GetVertNodeMatches: No input normal!"
        << "\n\t Returning empty result.";
        return matches;
    }

    double align_tol = 0.0;
    double min_abs_dist = DBL_MAX;
    double min_aligned_dist = DBL_MAX;
    size_t closest_abs_node_idx = -1;
    size_t closest_node_idx = -1;

    // find nearest node to point
    for (size_t q_pt = 0; q_pt < query_cloud->size(); q_pt++)
    {
        double sq_dist = 0.0;
        sq_dist += (source.x - query_cloud->points[q_pt].x)*(source.x - query_cloud->points[q_pt].x);
        sq_dist += (source.y - query_cloud->points[q_pt].y)*(source.y - query_cloud->points[q_pt].y);
        sq_dist += (source.z - query_cloud->points[q_pt].z)*(source.z - query_cloud->points[q_pt].z);
        double dist = sqrt(sq_dist);

        Eigen::Vector3d q_n{
            query_cloud->points[q_pt].normal_x,
            query_cloud->points[q_pt].normal_y,
            query_cloud->points[q_pt].normal_z
        };

        if (utils::areAligned(src_n,q_n,align_tol) && dist <= min_aligned_dist)
        {
            min_aligned_dist = dist; 
            closest_node_idx = q_pt;
        }

        //plan b
        if (dist < min_abs_dist)
        {
            min_abs_dist = dist;
            closest_abs_node_idx = q_pt;
        }

    }
    
    if (closest_node_idx == -1)
    {
        std::string e 
            = "DeformGraph::GetVertNodeMatches: error finding nearest node";
        BOOST_LOG_TRIVIAL(error) << e;
        throw std::runtime_error(e);
    }

    // if nearest aligned point is too far just take nearest
    if (min_aligned_dist > this->node_sample_radius_)
    {
        closest_node_idx = closest_abs_node_idx;
    }
    
    matches.push_back(closest_node_idx);

    // add 1-ring connected neighborhood
    std::vector<int> one_ring = this->node_neigh[closest_node_idx];
    for (size_t nid = 0; nid < one_ring.size(); nid++)
    {
        matches.push_back(one_ring[nid]);
    }

    // add 2-ring connected neighborhood
    for (size_t nid = 0; nid < one_ring.size(); nid++)
    {
        std::vector<int> tmp_ring = this->node_neigh[one_ring[nid]];
        for (size_t tid = 0; tid < tmp_ring.size(); tid++)
        {
            matches.push_back(tmp_ring[tid]);
        }
    }

    // remove duplicates (unsorted!)
    matches.erase( 
        utils::remove_duplicates( 
            matches.begin(), matches.end()), 
        matches.end());

    // init vector pair: dist,idx
    typedef std::pair<double,int> npair;
    std::vector<npair> dist_vid;
    for (size_t match_id = 0; match_id < matches.size(); match_id++)
    {
        int mid = matches[match_id];
        double sq_dist = 0.0;
        sq_dist += (source.x - query_cloud->points[mid].x)*(source.x - query_cloud->points[mid].x);
        sq_dist += (source.y - query_cloud->points[mid].y)*(source.y - query_cloud->points[mid].y);
        sq_dist += (source.z - query_cloud->points[mid].z)*(source.z - query_cloud->points[mid].z);
        double dist  = sqrt(sq_dist);
        
        if(dist < (1.5*this->node_sample_radius_*this->dmax_scaler_))
        {
            dist_vid.push_back(std::make_pair(dist,matches[match_id]));
        }
    }

    // sort by distance
    std::sort(dist_vid.begin(), dist_vid.end(),
    [](const npair& l, const npair& r) {
        if (l.first != r.first){
        return l.first < r.first;}
        return l.second < r.second; //unlikely to need
    });   

    std::vector<int> output;
    for (size_t match_id = 0; match_id < dist_vid.size(); match_id++)
    {
        output.push_back(dist_vid[match_id].second);
    }

    if (output.empty())
    {
        // Throw error and default to global
        std::string e 
            = "DeformGraph::GetVertNodeMatches: No graph node met criteria for point";
        BOOST_LOG_TRIVIAL(error) << e;
        throw std::runtime_error(e);

    }
    

    return output;
}
////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<int>> DeformGraph::GetConNodeNeighbours(
    const std::vector<pcl::PointXYZRGBNormal>& cons,
    int search_size)
{
    std::vector<std::vector<int>> node_knn_lists;

    Matcher mt;
    MeshProcessing mp;

    CloudPtr query_cloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    pcl::fromPCLPointCloud2(this->graph_structure_.cloud,*query_cloud);

    for (size_t cid = 0; cid < cons.size(); cid++)
    {
        int search_size = kDG_KNN;
        std::vector<int> nodes;

        double min_dist = DBL_MAX;
        size_t closest_node_idx = 0;

        // brute force, EVERY point in destination cloud
        // find nearest node to points
        for (size_t q_pt = 0; q_pt < query_cloud->size(); q_pt++)
        {
            double sq_dist = 0.0;
            sq_dist += (cons[cid].x - query_cloud->points[q_pt].x)*(cons[cid].x - query_cloud->points[q_pt].x);
            sq_dist += (cons[cid].y - query_cloud->points[q_pt].y)*(cons[cid].y - query_cloud->points[q_pt].y);
            sq_dist += (cons[cid].z - query_cloud->points[q_pt].z)*(cons[cid].z - query_cloud->points[q_pt].z);
            double dist = sqrt(sq_dist);

            if ( dist <= min_dist)
            {
            min_dist = dist; 
            closest_node_idx = q_pt;
            }
        }

        // find neighbours by 1-ring connectivity
        nodes = this->node_neigh[closest_node_idx];

        // sort by distance.....
        typedef std::pair<double,int> npair;
        std::vector<npair> dist_vid;
        for (size_t node_id = 0; node_id < nodes.size(); node_id++)
        {
            int mid = nodes[node_id];
            double sq_dist = 0.0;
            sq_dist += (cons[cid].x - query_cloud->points[mid].x)*(cons[cid].x - query_cloud->points[mid].x);
            sq_dist += (cons[cid].y - query_cloud->points[mid].y)*(cons[cid].y - query_cloud->points[mid].y);
            sq_dist += (cons[cid].z - query_cloud->points[mid].z)*(cons[cid].z - query_cloud->points[mid].z);
            double dist  = sqrt(sq_dist);
            dist_vid.push_back(std::make_pair(dist,nodes[node_id]));
        }

        std::sort(dist_vid.begin(), dist_vid.end(),
        [](const npair& l, const npair& r) {
            if (l.first != r.first){
            return l.first < r.first;}
            return l.second < r.second; //unlikely to need
        });   

        std::vector<int> sorted_nodes;
        for (size_t node_id = 0; node_id < dist_vid.size(); node_id++)
        {
            sorted_nodes.push_back(dist_vid[node_id].second);
        }

        node_knn_lists.push_back(sorted_nodes);
    }

    return node_knn_lists;
}

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::NormalizeTforms()
{

    for (size_t n = 0; n < this->node_trans.size(); n++)
    {
        Eigen::Vector3d summed_tform; 
        summed_tform.setZero();
        int num_neigh = std::ceil(kDG_KNN); 
        for (size_t nn = 0; nn < num_neigh; nn++)
        {
            summed_tform += node_trans[node_neigh[n][nn]];
        }
        summed_tform /= num_neigh;
        node_trans[n] = summed_tform;
    }
}

////////////////////////////////////////////////////////////////////////////////
void DeformGraph::Deform(
    const std::vector<std::vector<int>>& node_influence_list,
    const std::vector<std::vector<double>>& vert_node_weights,
    const pcl::PolygonMesh& input_graph,
    TriMesh& d_mesh) 
{

    // get pt and norm info
    MeshProcessing mp;
    CloudPtr graph_cloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    pcl::PolygonMesh tmp_mesh = input_graph;
    graph_cloud = mp.GetCloudFromPolygonMesh(tmp_mesh);

	std::vector<double> weights;
    std::vector<Eigen::Vector3d>& verts = d_mesh.vertex_coord;
    std::vector<Eigen::Vector3d>& norms = d_mesh.norm_coord;
    std::vector<Eigen::Vector3d> vs; 
    
	for (size_t i=0; i<verts.size(); ++i) 
    {

        std::vector<double> weights = vert_node_weights[i];

        Eigen::Vector3d vert = verts[i], norm = norms[i]; 
		Eigen::Matrix3d rot = Eigen::Matrix3d::Zero(); 
        verts[i] = Eigen::Vector3d(0.0, 0.0, 0.0);

        for (int j=0; j<node_influence_list[i].size(); ++j) 
        {
            int cur_node_idx = node_influence_list[i][j];
            Eigen::Vector3d g_j = this->node_pos[cur_node_idx];

            Eigen::Matrix3d R_j = this->node_rots_mat[cur_node_idx];
            Eigen::Vector3d t_j = this->node_trans[cur_node_idx];
            
            verts[i] += weights[j] * (R_j*(vert - g_j) + g_j + t_j);  
        
        }    
    }
}

////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd DeformGraph::InitAffineVector(size_t num_nodes) const
{
    Eigen::VectorXd affine_vector(12 * num_nodes);

    #pragma omp parallel for
    for (int node_idx = 0; node_idx < num_nodes; ++node_idx)
    {
        // R
        affine_vector(12 * node_idx + 0) = 1.0;
        affine_vector(12 * node_idx + 1) = 0.0;
        affine_vector(12 * node_idx + 2) = 0.0;
        
        affine_vector(12 * node_idx + 3) = 0.0;
        affine_vector(12 * node_idx + 4) = 1.0;
        affine_vector(12 * node_idx + 5) = 0.0;

        affine_vector(12 * node_idx + 6) = 0.0;
        affine_vector(12 * node_idx + 7) = 0.0;
        affine_vector(12 * node_idx + 8) = 1.0;

        // t
        affine_vector(12 * node_idx + 9) = 0.0;
        affine_vector(12 * node_idx + 10) = 0.0;
        affine_vector(12 * node_idx + 11) = 0.0;
    }
    
    return affine_vector;
}

////////////////////////////////////////////////////////////////////////////////

double DeformGraph::OptimizeOnce(
    DeformGraph& graph,   
    GNParams& params,
    std::vector<std::vector<int>>& con_node_neigh,
    std::vector<std::vector<double>>& con_node_weights,
	const NodeCorrs& cons) 
{

    // Build Stacked Affine Vector x; 
    Eigen::VectorXd x = InitAffineVector(graph.node_pos.size());

    GaussNewtonSolver gn_solver;
    gn_solver.SetParams(params);

    double e_min = gn_solver.solve(
        graph.node_neigh,
        graph.node_node_weights,
        cons,
        con_node_neigh,
        con_node_weights,
        graph.node_pos,
        x);

    // update graph
    for (size_t node_idx = 0; node_idx < graph.node_pos.size(); node_idx++)
    {
        int vn = node_pos.size();
        std::vector<Eigen::Matrix3d> rots(vn); 
        std::vector<Eigen::Vector3d> trans(vn);

        for (int i=0; i<vn; ++i) 
        {
            for (int r=0; r<3; ++r) 
            {
                for (int c=0; c<3; ++c)
                {
                    rots[i](r, c) = x[12*i+3*r+c];
                }
                trans[i][r] = x[12*i+9+r];
            }
        }
        for (int i=0; i<vn; ++i) 
        {
            node_rots_mat[i] = rots[i];
            Eigen::Quaterniond q_rot(rots[i]);
            node_rots[i] = q_rot;
            node_trans[i] = trans[i];
        }

    }
    
	return e_min;
}

////////////////////////////////////////////////////////////////////////////////

void DeformGraph::DebugRandomVertNodeDeformInfluence(
  std::vector<std::vector<int>> node_influence_list,
  std::vector<std::vector<double>> vert_node_weights,
  CloudPtr moving_cloud,
  pcl::PolygonMesh graph_mesh,
  std::string out_path
)
{

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, node_influence_list.size());

    int random = debug_vertex_;//dis(gen);

    // init all gray
    MeshProcessing mp;
    CloudPtr input_graph_cloud = mp.GetCloudFromPolygonMesh(graph_mesh);
    for (size_t pt = 0; pt < input_graph_cloud->size(); pt++)
    {
    input_graph_cloud->points[pt].r=50;
    input_graph_cloud->points[pt].g=50;
    input_graph_cloud->points[pt].b=50;
    }

    // Get K ring neighbours
    std::vector<int> nn_nodes = node_influence_list[random];
    std::cerr<<"debug index: "<<random<<"\n";
    std::cerr<<"neighbours: "<<nn_nodes.size()<<'\n';
    std::vector<double> weights = vert_node_weights[random];

    auto max_it = std::max_element(weights.begin(),weights.end());
    double max_w(*max_it);

    // init neighbors green
    for (size_t nn = 0; nn < nn_nodes.size(); nn++)
    {
        double tmp_w = weights[nn]/max_w;
        double color = utils::LERP(0,255,tmp_w);
        // std::cerr<<"["<<nn_nodes[nn]<<": "<<color<<"]";
        input_graph_cloud->points[nn_nodes[nn]].r=0;
        input_graph_cloud->points[nn_nodes[nn]].g=color;
        input_graph_cloud->points[nn_nodes[nn]].b=0;
    }
    std::cerr<<"\n";

    // // init seed blue
    // input_graph_cloud->points[nn_nodes[0]].r=0;
    // input_graph_cloud->points[nn_nodes[0]].g=0;
    // input_graph_cloud->points[nn_nodes[0]].b=255;

    // Add selected point in red 
    pcl::PointXYZRGBNormal pt(moving_cloud->points[random]);
    pt.r = 255;
    pt.g = 0;
    pt.b = 0;
    input_graph_cloud->push_back(pt);

    pcl::io::savePLYFile(out_path,*input_graph_cloud);
    std::cerr<<"finished kNN test"<<'\n';
    // std::cin.get();
    // exit(1);
}

////////////////////////////////////////////////////////////////////////////////