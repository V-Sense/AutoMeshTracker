
#include "mesh_processing.h"
 

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::CreateSmoothedAbstractionLayer(
    const pcl::PolygonMesh& in_mesh,
    const int& target_vertices,
    const int decim_type)
{

    pcl::PolygonMesh output;

    switch (decim_type)
    {
    case MeshProcDefs::PQEC:
      output = PoissonQECAbstraction(in_mesh,target_vertices);
      break;
    case MeshProcDefs::VG:
    // TODO: Stub only, needs implementation
      BOOST_LOG_TRIVIAL(debug) 
        << "Warning, VoxelGrid abstraction not implemented yet!";
      output = VGAbstraction(in_mesh,target_vertices);
      break;
    }

    // remove small components and ensure manifoldness
    output = RemoveSmallComponents(output,0.01);
    output = FixManifoldness(output);

    // recalculate normals
    output = CalculateVertexNormals(output);

    return output;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::PoissonQECAbstraction(
  const pcl::PolygonMesh& in_mesh,
  const int& target_vertices
)
{
    // BOOST_LOG_TRIVIAL(debug) << "Creating abstraction layer via QEC";
    pcl::PolygonMesh output = in_mesh;

  // clean input
  output = RemoveSmallComponents(output,0.01);
  output = FixManifoldness(output);

  // // decimate result
  output = DecimateMesh(output,target_vertices);
  double target_edge_length = 0.5*CalculateAverageEdgeLength(output);
  output = IsotropicRemesh(output,target_edge_length);

  // Laplacian smoothing
  pcl::MeshSmoothingLaplacianVTK vtk;
  boost::shared_ptr<pcl::PolygonMesh> in_mesh_vtk(new pcl::PolygonMesh(output));
  vtk.setInputMesh(in_mesh_vtk);
  vtk.setNumIter(100); 
  vtk.setFeatureEdgeSmoothing(true); 
  vtk.setFeatureAngle(M_PI/2);
  vtk.process(output);
  return output;

}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::VGAbstraction(
  const pcl::PolygonMesh& in_mesh,
  const int& target_vertices)
{

  pcl::PolygonMesh output(in_mesh);
  // TODO: implementation
  return output;

}


////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::SmoothMesh(
  const pcl::PolygonMesh& in_mesh,
  const int& iters,
  const bool edge_smoothing,
  const float feat_angle
)
{

    // BOOST_LOG_TRIVIAL(debug) << "smoothing...";
    pcl::PolygonMesh output = in_mesh;
    pcl::MeshSmoothingLaplacianVTK vtk;
    boost::shared_ptr<pcl::PolygonMesh> in_mesh_vtk(new pcl::PolygonMesh(output));
    vtk.setInputMesh(in_mesh_vtk);
    vtk.setNumIter(iters); 
    // vtk.setFeatureEdgeSmoothing(true); 
    // vtk.setFeatureAngle(M_PI/2);
    vtk.process(output);
    
    return output;
}


////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::DecimateMesh(
	const pcl::PolygonMesh& in_mesh,
	const int& target_verts
)
{
	pcl::PolygonMesh output = in_mesh;
	// BOOST_LOG_TRIVIAL(debug) << "decimating...";
	double in_size = (double)output.cloud.width;
	const double r_factor = 1.0-(target_verts/in_size);
	const double pct = r_factor*100;
	// BOOST_LOG_TRIVIAL(debug) << "decimating to "+std::to_string(pct)+"%...";

	pcl::MeshQuadricDecimationVTK vtk;
	boost::shared_ptr<pcl::PolygonMesh> in_mesh_vtk(new pcl::PolygonMesh(output));
	vtk.setInputMesh(in_mesh_vtk);
	vtk.setTargetReductionFactor(r_factor);
	vtk.process(output);

	output = FixManifoldness(output);
  return output;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::PoissonRemesh(
    const pcl::PolygonMesh& in_mesh,
    const int& oct_depth
)
{
    // get cloud data
    pcl::PolygonMesh tmp(in_mesh);
    CloudPtr mesh_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
    mesh_cloud = GetCloudFromPolygonMesh(tmp);

    pcl::Poisson<pcl::PointXYZRGBNormal> poisson;
    poisson.setDepth(oct_depth);
    poisson.setInputCloud(mesh_cloud);

    pcl::PolygonMesh outmesh;
    poisson.reconstruct(outmesh);
    return outmesh;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::IsotropicRemesh(
  const pcl::PolygonMesh& in_mesh,
  double target_edge_length)
  {
    
    // default condition
    if (target_edge_length == 0)
    {
      // calculate the average edge length of input mesh
      target_edge_length = 0.75 * CalculateAverageEdgeLength(in_mesh);
    }

    std::stringstream iss = utils::PclMesh2CgalStream(in_mesh);

    Surface_Mesh_3 mesh;
    if (!iss || !(iss >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
      std::string e = "In MeshProcessing::IsotropicRemesh: invalid input";
      BOOST_LOG_TRIVIAL(fatal) << e;
      throw std::runtime_error(e);
    }

    unsigned int nb_iter = 3;
      std::vector<sm3_edge_descriptor> border;
      PMP::border_halfedges(faces(mesh),
        mesh,
        boost::make_function_output_iterator(halfedge2edge(mesh, border)));
      PMP::split_long_edges(border, target_edge_length, mesh);
    PMP::isotropic_remeshing(
        faces(mesh),
        target_edge_length,
        mesh,
        PMP::parameters::number_of_iterations(nb_iter)
        .protect_constraints(true)
        );
    
    std::stringstream oss; 
    oss << mesh;
    pcl::PolygonMesh out_mesh = utils::CgalStream2PclMesh(oss);
    return out_mesh;
  }
  
////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::RemeshInstantMeshes(
  const pcl::PolygonMesh& in_mesh,
  const size_t& target_vertices)
{

  pcl::PolygonMesh outmesh;
  outmesh = PoissonRemesh(in_mesh,9);
  outmesh = FixManifoldness(outmesh);
  outmesh = RemoveSmallComponents(outmesh,0.01);

  // save
  std::string f_tmp_out = "/tmp/in_mesh.ply";
  pcl::io::savePLYFile(f_tmp_out,outmesh);

  // meshlab poisson and remesh
  std::string f_remeshed = "/tmp/remeshed.ply";
  std::string remesh_cmd 
    = "/MeshTracker/resource/thirdparty/Instant_Meshes -o "+
    f_remeshed+" -d -i -v "+std::to_string(target_vertices)+" -r 6 -p 6 "+
    f_tmp_out;

  system(remesh_cmd.c_str());
  pcl::io::loadPolygonFilePLY(f_remeshed,outmesh);

  outmesh = FixManifoldness(outmesh);
  outmesh = RemoveSmallComponents(outmesh,0.01);

  return outmesh;

}

////////////////////////////////////////////////////////////////////////////////

// pcl::PolygonMesh MeshProcessing::PoissonMeshing(
//     const CloudPtr& in_cloud,
//     const int& oct_depth
// )
// {
//     // Get normals
//     CloudPtr cloud_norm (new pcl::PointCloud<pcl::PointXYZRGBNormal>(*in_cloud));
//     CloudProcessing cp;
//     cp.CalculatePointCloudNormals(cloud_norm);

//     pcl::Poisson<pcl::PointXYZRGBNormal> poisson;
//     poisson.setDepth(oct_depth);
//     poisson.setInputCloud(cloud_norm);

//     pcl::PolygonMesh out_mesh;
//     poisson.reconstruct(out_mesh);
//     return out_mesh;
// }

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::RemoveSmallComponents
(
    const pcl::PolygonMesh& in_mesh, 
    const double& percent
)
{
    Mesh he_mesh_in;
    pcl::geometry::toHalfEdgeMesh(in_mesh,he_mesh_in);

    std::vector<bool> visited(he_mesh_in.sizeVertices(), false);
    std::vector<std::vector<Mesh::VertexIndex>> components;

    unsigned int nVtx = he_mesh_in.sizeVertices();
    unsigned int nComponents = 0;
    Mesh::VertexIndex vi;

    while (true)
    {
        // Find an unvisited vertex
        bool found = false;
        for (unsigned int i = 0; i < he_mesh_in.sizeVertices(); i++)
        {
            if(!visited[i])
            {
                found = true;
                visited[i] = true;
                vi = Mesh::VertexIndex(i);
                break;
            }
        }

        // If none was found -> finished
        if (!found)
        {
            break;
        } 
        nComponents++;

        std::vector<Mesh::VertexIndex> vindices;
        std::vector<Mesh::VertexIndex> componentIndices;
        vindices.push_back(vi);

        while(vindices.size() > 0)
        {
            // check the end of the stack and pop it out
            Mesh::VertexIndex current = vindices.back();
            if (!current.isValid())
            {
                continue;
            } 
            vindices.pop_back();

            Mesh::VertexAroundVertexCirculator vac = he_mesh_in.getVertexAroundVertexCirculator(current);
            const Mesh::VertexAroundVertexCirculator vac_end = vac;
            if(!vac.isValid())
            {
                continue;
            }

            // check neighboring vertices and see if they have been already visited
            do 
            {
                if (!visited[vac.getTargetIndex().get()])
                {
                    visited[vac.getTargetIndex().get()] = true;
                    vindices.push_back(vac.getTargetIndex());
                    componentIndices.push_back(vac.getTargetIndex());
                }

            } while (++vac != vac_end);
        }

        components.push_back(componentIndices);

    }

    if (components.size() != 1)
    {
        for (unsigned int i = 0; i < components.size(); i++)
        {
            if (components[i].size() > (int) (nVtx * percent))
            {
                continue;
            }

            const std::vector<Mesh::VertexIndex> curr_comp = components[i];
            for (unsigned int j = 0; j < curr_comp.size(); j++)
            {
                const Mesh::VertexIndex curr_vi = curr_comp[j];
                he_mesh_in.deleteVertex(curr_vi);
            }
        }
    }

    he_mesh_in.cleanUp();
    pcl::PolygonMesh out_mesh;
    pcl::geometry::toFaceVertexMesh(he_mesh_in,out_mesh);
    return out_mesh;
}

////////////////////////////////////////////////////////////////////////////////

double MeshProcessing::CalculateAverageEdgeLength(const pcl::PolygonMesh& mesh)
{

  Eigen::MatrixXd V;
	Eigen::MatrixXi F;

  typedef std::pair<int,int> edge;
  std::vector<edge> edges; 

	Mesh he_mesh;
  pcl::geometry::toHalfEdgeMesh(mesh, he_mesh);

	size_t num_verts = he_mesh.sizeVertices();
	size_t num_faces = he_mesh.sizeFaces();

	V.setConstant(num_verts,3,0);
	F.setConstant(num_faces,3,0); // trimesh only!
	Mesh::VertexDataCloud v_cloud(he_mesh.getVertexDataCloud());

	for (size_t vert = 0; vert < num_verts; vert++)
	{	
		V(vert,0) = v_cloud.points[vert].x;
		V(vert,1) = v_cloud.points[vert].y;
		V(vert,2) = v_cloud.points[vert].z;
	}

  // Assuming trimesh
	for (size_t face = 0; face < num_faces; face++)
	{
		Mesh::FaceIndex fi(face);
		Mesh::VertexAroundFaceCirculator vafc;
		vafc = he_mesh.getVertexAroundFaceCirculator(fi);
		if(!vafc.isValid())
		{
			continue;
		}

    int v0 = vafc.getTargetIndex().get(); ++vafc;
    int v1 = vafc.getTargetIndex().get(); ++vafc;
    int v2 = vafc.getTargetIndex().get();

    // always store edges as smallest index first
    edge e1 = v0 < v1 ? std::make_pair(v0,v1) : std::make_pair(v1,v0);
    edge e2 = v0 < v2 ? std::make_pair(v0,v2) : std::make_pair(v2,v0);
    edge e3 = v1 < v2 ? std::make_pair(v1,v2) : std::make_pair(v2,v1);

    edges.push_back(e1);
    edges.push_back(e2);
    edges.push_back(e3);
	}

  // remove duplicates
  std::sort(edges.begin(),edges.end());
  edges.erase(std::unique(
      edges.begin(), edges.end(), utils::PairSortPredicate), edges.end());

  double total_len = 0.0;
  for (size_t edge = 0; edge < edges.size(); edge++)
  {
    total_len += std::abs(V(edges[edge].first,0) - V(edges[edge].second,0));
    total_len += std::abs(V(edges[edge].first,1) - V(edges[edge].second,1));
    total_len += std::abs(V(edges[edge].first,2) - V(edges[edge].second,2)); 
  }
  
  return total_len/edges.size();
}

////////////////////////////////////////////////////////////////////////////////

double MeshProcessing::GetAffineScaleFactor
(
    const pcl::PolygonMesh& in_mesh, 
    const double& target_height,
    const bool& translate2origin
)
{
    pcl::PolygonMesh out_mesh = in_mesh;
    CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::fromPCLPointCloud2(out_mesh.cloud,*cloud);

    // Get bounding box
    mesh_processing::BBox bb(out_mesh);

    // assuming height refers to longest axis, get height scale factor
    double height = std::max({
        bb.max_x - bb.min_x,
        bb.max_y - bb.min_y,
        bb.max_z - bb.min_z});
    
    return target_height/height;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::ApplyScaleFactor
(
  const pcl::PolygonMesh     & in_mesh, 
  const double               &      sf,
  const std::array<double,3> &  offset
)
{

    pcl::PolygonMesh out_mesh = in_mesh;
    CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::fromPCLPointCloud2(out_mesh.cloud,*cloud);

    // Get bounding box
    mesh_processing::BBox bb(out_mesh);


    // translate center to world origin and store
    for (size_t pt = 0; pt < cloud->points.size(); pt++)
    {
        cloud->points[pt].x -= offset[0];
        cloud->points[pt].y -= offset[1];
        cloud->points[pt].z -= offset[2];
    }

    // apply factor to each dimension 
    // translate center to world origin and store
    for (size_t pt = 0; pt < cloud->points.size(); pt++)
    {
        cloud->points[pt].x *= sf;
        cloud->points[pt].y *= sf;
        cloud->points[pt].z *= sf;
    }

    // invert translation 
    // for (size_t pt = 0; pt < cloud->points.size(); pt++)
    // {
    //     cloud->points[pt].x += offset[0];
    //     cloud->points[pt].y += offset[1];
    //     cloud->points[pt].z += offset[2];
    // }
    
    pcl::toPCLPointCloud2(*cloud,out_mesh.cloud);
    return out_mesh;

}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::FixManifoldness(const pcl::PolygonMesh &mesh)
{

  // Convert to half edge mesh
  Mesh he_mesh;
  pcl::geometry::toHalfEdgeMesh(mesh, he_mesh);

  for (unsigned int i = 0; i < he_mesh.sizeVertices(); i++)
  {
    Mesh::VertexIndex vi = Mesh::VertexIndex(i);
    if (!vi.isValid())
    {
      he_mesh.deleteVertex(vi);
      continue;
    }

    if (!he_mesh.isManifold(vi))
    {

      Mesh::FaceAroundVertexCirculator fav = he_mesh.getFaceAroundVertexCirculator(vi);
      if (!fav.isValid())
        continue;
      const Mesh::FaceAroundVertexCirculator fav_end = fav;
      std::vector<Mesh::FaceIndex> toDelete;
      do
      {
        Mesh::FaceIndex fi = fav.getTargetIndex();
        toDelete.push_back(fi);
      } while (++fav != fav_end);

      for (unsigned int j = 0; j < toDelete.size(); j++)
      {
        if (toDelete[j].isValid())
        {
          he_mesh.deleteFace(toDelete[j]);
        }
      }
    }
  }

  he_mesh.cleanUp();
  pcl::PolygonMesh outmesh;

  pcl::geometry::toFaceVertexMesh(he_mesh, outmesh);
  return outmesh;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PointXYZRGBNormal MeshProcessing::GetOneRingCentroid(
    const pcl::PolygonMesh &mesh,
    const int &init_pt)
{

  // convert to geometry processing format
  Mesh he_mesh;
  pcl::geometry::toHalfEdgeMesh(mesh, he_mesh);

  // extract point cloud
  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr mesh_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(mesh.cloud, *mesh_cloud);

  // initialise output
  pcl::PointXYZRGBNormal centroid = mesh_cloud->points[init_pt];

  // initialize the vert-around-vert circulator
  Mesh::VertexIndex vidx = Mesh::VertexIndex(init_pt);
  Mesh::VertexAroundVertexCirculator vavc = he_mesh.getVertexAroundVertexCirculator(vidx);

  // define endpoint
  const Mesh::VertexAroundVertexCirculator vavc_end = vavc;

  // will throw errors otherwise
  if (!vavc.isValid())
  {
    BOOST_LOG_TRIVIAL(warning) << "In GetOneRingCentroid: vav circulator not valid. \nReturning original point ";
    return centroid;
  };

  // sum each vertex around center point and normalize
  Eigen::Vector3f addedpoints(0.0, 0.0, 0.0);
  int pts = 0;
  do
  {

    int idx = vavc.getTargetIndex().get();
    addedpoints(0) += mesh_cloud->points[idx].x;
    addedpoints(1) += mesh_cloud->points[idx].y;
    addedpoints(2) += mesh_cloud->points[idx].z;

    pts++;

  } while (++vavc != vavc_end);

  // assign new data
  centroid.x = addedpoints(0) / pts;
  centroid.y = addedpoints(1) / pts;
  centroid.z = addedpoints(2) / pts;

  return centroid;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshProcessing::GetOneRingNeighbors(
    const pcl::PolygonMesh &mesh,
    const int &init_pt)
{

  // convert to geometry processing format
  Mesh he_mesh;
  pcl::geometry::toHalfEdgeMesh(mesh, he_mesh);

  // extract point cloud
  pcl::PointCloud<pcl::PointXYZRGBNormal>::Ptr mesh_cloud
    (new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(mesh.cloud, *mesh_cloud);

  // initialise output
  std::vector<int> neighbours;

  // initialize the vert-around-vert circulator
  Mesh::VertexIndex vidx = Mesh::VertexIndex(init_pt);
  Mesh::VertexAroundVertexCirculator vavc = 
    he_mesh.getVertexAroundVertexCirculator(vidx);

  // define endpoint
  const Mesh::VertexAroundVertexCirculator vavc_end = vavc;

  // will throw errors otherwise
  if (!vavc.isValid())
  {
    BOOST_LOG_TRIVIAL(warning)
      <<"In GetOneRingNeighbors: vav circulator not valid.\n"
      <<"Returning empty list ";
    return neighbours;
  };

  do
  {
    int idx = vavc.getTargetIndex().get();
    if (idx == init_pt)
    {
      continue;
    }
    
    neighbours.push_back(idx);
  } 
  while(++vavc != vavc_end);

  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> MeshProcessing::GetTwoRingNeighbors(
  const pcl::PolygonMesh& mesh,
  const int& vid)
{
  

  // get first ring
  std::vector<int> neighbours = GetOneRingNeighbors(mesh,vid);
  size_t last_size = neighbours.size();
  std::vector<int> sizes;
  
  // get next ring neighbours
  int start_idx = 0;


  std::vector<int> nn;
  int n = start_idx;
  do
  {
    nn.clear();
    nn = GetOneRingNeighbors(mesh,neighbours[n]);
    
    // add new neighbours
    neighbours.insert( neighbours.end(), nn.begin(), nn.end() );
    ++n;
  }
  while(n < last_size);

  // remove duplicates (unsorted!)
  neighbours.erase( 
    utils::remove_duplicates( 
      neighbours.begin(), 
      neighbours.end() ), 
      neighbours.end() );

  last_size = neighbours.size();
  sizes.push_back(last_size);

  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////


std::vector<int> MeshProcessing::GetKRingNeighbors(
  const pcl::PolygonMesh& mesh,
  const int& vid,
  const int& k)
{
  if (k == 1 )
  {
    return GetOneRingNeighbors(mesh,vid);
  }

  if (k < 2)
  {
    std::string e = "In GetKRingNeighbours: invalid k value -> "+k;
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  // size_t max_nn = mesh.cloud.width;
  // TODO: catch for large k

  // get first ring
  std::vector<int> neighbours = GetOneRingNeighbors(mesh,vid);
  size_t last_size = neighbours.size();
  std::vector<int> sizes;
  
  // get k ring neighbours
  int start_idx = 0;
  for (size_t r = 1; r < k; r++)
  {

    // set a skip index every odd iteration as everything from 2 iters prior
    // is already guaranteed 
    if(r%2==1 && r>2)
    {
      start_idx = sizes[r-3];
    }

    std::vector<int> nn;
    int n = start_idx;
    do
    {
      nn.clear();
      nn = GetOneRingNeighbors(mesh,neighbours[n]);
      
      // add new neighbours
      neighbours.insert( neighbours.end(), nn.begin(), nn.end() );
      ++n;
    }
    while(n < last_size);

    // remove duplicates (unsorted!)
    neighbours.erase( 
      utils::remove_duplicates( 
        neighbours.begin(), 
        neighbours.end() ), 
        neighbours.end() );

    last_size = neighbours.size();
    sizes.push_back(last_size);
  }
  return neighbours;
}


////////////////////////////////////////////////////////////////////////////////

bool MeshProcessing::NormalsExist(
    const pcl::PolygonMesh &mesh)
{

  // copy data from input
  CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  for (size_t pt = 0; pt < cloud->size(); pt++)
  {
    if (
      cloud->points[pt].normal_x == 0 &&
      cloud->points[pt].normal_y == 0 &&
      cloud->points[pt].normal_z == 0)
    {
      return false;
    }
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::CalculateVertexNormals(
  const  pcl::PolygonMesh &mesh_)
{

  std::stringstream ss_in = utils::PclMesh2CgalStream(mesh_);
  Surface_mesh mesh;
  ss_in >> mesh;
  
  if (mesh.is_empty())
  {
    std::string e = 
      "In MeshProcessing::CalculateVertexNormals, couldn't read input file";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }

  auto fnormals = mesh.add_property_map<face_descriptor, Vector>
      ("f:normals", CGAL::NULL_VECTOR).first;
  auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>
      ("v:normals", CGAL::NULL_VECTOR).first;
  CGAL::Polygon_mesh_processing::compute_normals(
        mesh,
        vnormals,
        fnormals,
        CGAL::Polygon_mesh_processing::parameters::vertex_point_map(mesh.points()).
        geom_traits(K()));

  pcl::PolygonMesh outmesh(mesh_);
  CloudPtr tmp_cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(mesh_.cloud,*tmp_cloud);
  
  unsigned int idx = 0;
  for(vertex_descriptor vd: vertices(mesh)){
    tmp_cloud->points[idx].normal_x = vnormals[vd][0];
    tmp_cloud->points[idx].normal_y = vnormals[vd][1];
    tmp_cloud->points[idx].normal_z = vnormals[vd][2];
    idx++;
  }

  pcl::toPCLPointCloud2(*tmp_cloud,outmesh.cloud);
  
  return outmesh;
}

////////////////////////////////////////////////////////////////////////////////

CloudPtr MeshProcessing::GetCloudFromPolygonMesh(pcl::PolygonMesh &mesh)
{

  MeshProcessing mp;
  pcl::PolygonMesh tmp;
  if (!mp.NormalsExist(mesh))
  {
    tmp = mp.CalculateVertexNormals(mesh);
    CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::fromPCLPointCloud2(tmp.cloud, *cloud);

    return cloud;
  }
  else
  {
    // copy data from input
    CloudPtr cloud(new pcl::PointCloud<pcl::PointXYZRGBNormal>);
    pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

    return cloud;
  }
}

////////////////////////////////////////////////////////////////////////////////

bool MeshProcessing::HasNanData(const pcl::PolygonMesh& mesh)
{

  CloudPtr cloud (new pcl::PointCloud<pcl::PointXYZRGBNormal>());
  pcl::fromPCLPointCloud2(mesh.cloud, *cloud);

  for (size_t pt = 0; pt < cloud->size(); pt++)
  {
    if (
      std::isnan(cloud->points[pt].x) ||
      std::isnan(cloud->points[pt].y) ||
      std::isnan(cloud->points[pt].z))
    {
      return true;
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

pcl::PolygonMesh MeshProcessing::LERPMesh(
    const pcl::PolygonMesh source,
    const pcl::PolygonMesh target,
    const float step)
{

  // check that inputs have same num vertices
  if (source.cloud.width != target.cloud.width)
  {
    std::string e = "In MeshProcessing::LERPMesh: Inputs don't match";
    BOOST_LOG_TRIVIAL(error) << e;
    throw std::runtime_error(e);
  }
  
  // TODO: add topology check too

  // extract clouds
  pcl::PolygonMesh tmp_s = source;
  pcl::PolygonMesh tmp_t = target;
  CloudPtr s_cloud = GetCloudFromPolygonMesh(tmp_s);
  CloudPtr t_cloud = GetCloudFromPolygonMesh(tmp_t);

  // perform LERP
  CloudPtr r_cloud(s_cloud);
  for (size_t pt = 0; pt < s_cloud->size(); pt++)
  {
    r_cloud->points[pt].x += (t_cloud->points[pt].x - s_cloud->points[pt].x)*step;
    r_cloud->points[pt].y += (t_cloud->points[pt].y - s_cloud->points[pt].y)*step;
    r_cloud->points[pt].z += (t_cloud->points[pt].z - s_cloud->points[pt].z)*step;

    r_cloud->points[pt].normal_x += (t_cloud->points[pt].normal_x - s_cloud->points[pt].normal_x)*step;
    r_cloud->points[pt].normal_y += (t_cloud->points[pt].normal_y - s_cloud->points[pt].normal_y)*step;
    r_cloud->points[pt].normal_z += (t_cloud->points[pt].normal_z - s_cloud->points[pt].normal_z)*step;
  }

  
  // attach to output and return
  pcl::PolygonMesh output = source;
  pcl::toPCLPointCloud2(*r_cloud, output.cloud);
  return output;
}  

////////////////////////////////////////////////////////////////////////////////

// std::vector<int> MeshProcessing::SortByGeodesicDistance(
//   const pcl::PolygonMesh& mesh,
//   const int& src_pt_idx,
//   const std::vector<int>& pts_indxs)
// {
//   // 2d vector of indices and distances. Distance will be int as it's approximated
//   // in descrete vertex intervals
//   std::vector<std::vector<int>> idx_dist;
//   for (size_t pt = 0; pt < pts_indxs.size(); pt++)
//   {
//     idx_dist[pt][0] = pt;
//     idx_dist[pt][1] = GetGeodesicDistance(mesh,src_pt_idx,pts_indxs[pt]);
//   }

//   // TODO: new method to sort by 2nd column....
//   std::sort(idx_dist.begin(),idx_dist.end(),utils::colsort);
//   std::vector<int> out; 

//   for (size_t pt = 0; pt < pts_indxs.size(); pt++)
//   {
//     out.push_back(idx_dist[pt][0]);
//   }
  
//   return out;
// }

// ////////////////////////////////////////////////////////////////////////////////

// int MeshProcessing::GetGeodesicDistance(
//   const pcl::PolygonMesh& mesh,
//   const int& src_pt,
//   const int& dst_pt)
// {

//   Triangle_mesh tm;
//   utils::PclMesh2CgalStream(mesh) >> tm;

//   //property map for the distance values to the source set
//   Vertex_distance_map vertex_distance = 
//     tm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;
//   vertex_descriptor source = *(vertices(tm).first);
  
//   CGAL::Heat_method_3::estimate_geodesic_distances(tm, vertex_distance, source);
  
//   for(vertex_descriptor vd : vertices(tm)){
//     std::cout << vd << " ("<< tm.point(vd) << ")"
//               <<  " is at distance " << get(vertex_distance, vd) << std::endl;
//   }

//   int distance =0;
//   return distance;
// }

////////////////////////////////////////////////////////////////////////////////