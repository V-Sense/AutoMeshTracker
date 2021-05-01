#include "gauss_newton_solver.h"


/*
Based on Code from the following work:
Guo, Kaiwen, et al. "Robust non-rigid motion tracking and surface reconstruction using l0 regularization." Proceedings of the IEEE International Conference on Computer Vision. 2015.
Please cite if using this code for any research purposes. 

Received on April 6th, 2020. 
Please contact Dr. Lebin Yiu for original code
*/

////////////////////////////////////////////////////////////////////////////////

double GaussNewtonSolver::solve(
    const std::vector<std::vector<int>>& node_node_neigh,
    const std::vector<std::vector<double>>& node_node_weights,
    const NodeCorrs& constraints,
    const std::vector<std::vector<int>>& con_node_neigh,
    const std::vector<std::vector<double>>& con_node_weights,
    std::vector<Eigen::Vector3d>& node_pos,
    Eigen::VectorXd& affine_vector
)
{
    double gaussNewtonEnergy = 0.0;
    double gaussNewtonEnergyPrev = 0.0;
    size_t gaussNewtonIterCnt = 0;

    BOOST_LOG_TRIVIAL(debug) <<"Solver params: "
        <<params_.m_alpha_rigid<<", "
        <<params_.m_alpha_smooth<<", "
        <<params_.m_alpha_point<<", "
        <<params_.m_alpha_plane;

    while (true)
    {
        //------------------------------------------------------------------------
        //	Initialize numerical containers
        //------------------------------------------------------------------------
        std::vector<Eigen::Triplet<double>> JrList;
        std::vector<Eigen::Triplet<double>> AsList;
        std::vector<Eigen::Triplet<double>> ApList;
        std::vector<Eigen::Triplet<double>> AqList;

        Eigen::SparseMatrix<double> Jr;
        Eigen::SparseMatrix<double> As;
        Eigen::SparseMatrix<double> Ap;
        Eigen::SparseMatrix<double> Aq;

        Eigen::VectorXd fr;
        Eigen::VectorXd bs;
        Eigen::VectorXd bp;
        Eigen::VectorXd bq;

        Eigen::SparseMatrix<double> A0(12 * node_node_neigh.size(), 12 * node_node_neigh.size());
        Eigen::SparseMatrix<double> A1(12 * node_node_neigh.size(), 12 * node_node_neigh.size());
        Eigen::SparseMatrix<double> A2(12 * node_node_neigh.size(), 12 * node_node_neigh.size());
        Eigen::SparseMatrix<double> A3(12 * node_node_neigh.size(), 12 * node_node_neigh.size());

        Eigen::VectorXd b0(12 * node_node_neigh.size());
        Eigen::VectorXd b1(12 * node_node_neigh.size());
        Eigen::VectorXd b2(12 * node_node_neigh.size());
        Eigen::VectorXd b3(12 * node_node_neigh.size());

    #pragma omp sections
        {
            //	First independent section: RIGIDITY TERM
            //	Construct Jr, fr, A0, and b0
    #pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct fr
                //------------------------------------------------------------------------
                std::vector<double> frVector;
                for (int i = 0; i < node_node_neigh.size(); ++i)
                {
                    double rigid[6] = { 0.0 };
                    rigid[0] += affine_vector(12 * i + 0) * affine_vector(12 * i + 3);
                    rigid[0] += affine_vector(12 * i + 1) * affine_vector(12 * i + 4);
                    rigid[0] += affine_vector(12 * i + 2) * affine_vector(12 * i + 5);

                    rigid[1] += affine_vector(12 * i + 3) * affine_vector(12 * i + 6);
                    rigid[1] += affine_vector(12 * i + 4) * affine_vector(12 * i + 7);
                    rigid[1] += affine_vector(12 * i + 5) * affine_vector(12 * i + 8);

                    rigid[2] += affine_vector(12 * i + 0) * affine_vector(12 * i + 6);
                    rigid[2] += affine_vector(12 * i + 1) * affine_vector(12 * i + 7);
                    rigid[2] += affine_vector(12 * i + 2) * affine_vector(12 * i + 8);

                    rigid[3] += affine_vector(12 * i + 0) * affine_vector(12 * i + 0);
                    rigid[3] += affine_vector(12 * i + 1) * affine_vector(12 * i + 1);
                    rigid[3] += affine_vector(12 * i + 2) * affine_vector(12 * i + 2);
                    rigid[3] -= 1.0;

                    rigid[4] += affine_vector(12 * i + 3) * affine_vector(12 * i + 3);
                    rigid[4] += affine_vector(12 * i + 4) * affine_vector(12 * i + 4);
                    rigid[4] += affine_vector(12 * i + 5) * affine_vector(12 * i + 5);
                    rigid[4] -= 1.0;

                    rigid[5] += affine_vector(12 * i + 6) * affine_vector(12 * i + 6);
                    rigid[5] += affine_vector(12 * i + 7) * affine_vector(12 * i + 7);
                    rigid[5] += affine_vector(12 * i + 8) * affine_vector(12 * i + 8);
                    rigid[5] -= 1.0;

                    frVector.push_back(rigid[0]);
                    frVector.push_back(rigid[1]);
                    frVector.push_back(rigid[2]);
                    frVector.push_back(rigid[3]);
                    frVector.push_back(rigid[4]);
                    frVector.push_back(rigid[5]);
                }

                fr.resize(frVector.size());
                for (int i = 0; i < frVector.size(); ++i)
                {
                    fr(i) = frVector[i];
                }

                //------------------------------------------------------------------------
                //	Construct Jr
                //------------------------------------------------------------------------
                size_t jOffset = 0;
                for (int i = 0; i < node_node_neigh.size(); ++i)
                {
                    double x[9] = { 0.0 };
                    x[0] = affine_vector(12 * i + 0);
                    x[1] = affine_vector(12 * i + 1);
                    x[2] = affine_vector(12 * i + 2);
                    x[3] = affine_vector(12 * i + 3);
                    x[4] = affine_vector(12 * i + 4);
                    x[5] = affine_vector(12 * i + 5);
                    x[6] = affine_vector(12 * i + 6);
                    x[7] = affine_vector(12 * i + 7);
                    x[8] = affine_vector(12 * i + 8);

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[5]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, x[8]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[5]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, x[8]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 0, 2 * x[0]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 1, 2 * x[1]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 2, 2 * x[2]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 3, 2 * x[3]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 4, 2 * x[4]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 5, 2 * x[5]));
                    ++jOffset;

                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 6, 2 * x[6]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 7, 2 * x[7]));
                    JrList.push_back(Eigen::Triplet<double>(jOffset, 12 * i + 8, 2 * x[8]));
                    ++jOffset;
                }

                Jr.resize(jOffset, 12 * node_node_neigh.size());
                Jr.setFromTriplets(JrList.begin(), JrList.end());

                Eigen::SparseMatrix<double> JrTJr = Eigen::SparseMatrix<double>(Jr.transpose()) * Jr;
                A0 = params_.m_alpha_rigid * JrTJr;
                b0 = params_.m_alpha_rigid * JrTJr * affine_vector.cast<double>();
                b0 -= params_.m_alpha_rigid * Jr.transpose() * fr;
            }

            //	Second independent section: SMOOTHING TERM
            //	Construct As, bs, A1 and b1
    #pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct As and bs
                //------------------------------------------------------------------------
                std::vector<double> bsVector;
                size_t sOffset = 0;

                for (int node_idx = 0; node_idx < node_node_neigh.size(); ++node_idx)
                {
                    Eigen::Vector3d v0 = node_pos[node_idx];

                    for (int neigh_idx = 0; neigh_idx < node_node_neigh[node_idx].size(); ++neigh_idx)
                    {
                        double weight = node_node_weights[node_idx][neigh_idx];
                        // double weightRoot = sqrt(weight); //should be root already, keep jic
                        int nid_actual = node_node_neigh[node_idx][neigh_idx];
                        Eigen::Vector3d v1 = node_pos[nid_actual];

                        double vec[3] = { 0.0 };
                        vec[0] = weight * (v1[0] - v0[0]);
                        vec[1] = weight * (v1[1] - v0[1]);
                        vec[2] = weight * (v1[2] - v0[2]);

                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 0, vec[0]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 3, vec[1]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 6, vec[2]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 9, weight));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * nid_actual + 9, -weight));
                        ++sOffset;

                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 1, vec[0]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 4, vec[1]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 7, vec[2]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 10, weight));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * nid_actual + 10, -weight));
                        ++sOffset;

                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 2, vec[0]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 5, vec[1]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 8, vec[2]));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * node_idx + 11, weight));
                        AsList.push_back(Eigen::Triplet<double>(sOffset, 12 * nid_actual + 11, -weight));
                        ++sOffset;

                        bsVector.push_back(vec[0]);
                        bsVector.push_back(vec[1]);
                        bsVector.push_back(vec[2]);
                    }
                }

                As.resize(sOffset, 12 * node_node_neigh.size());
                As.setFromTriplets(AsList.begin(), AsList.end());

                bs.resize(bsVector.size());
                for (int i = 0; i < bsVector.size(); ++i)
                {
                    bs(i) = bsVector[i];
                }

                A1 = params_.m_alpha_smooth * Eigen::SparseMatrix<double>(As.transpose()) * As;
                b1 = params_.m_alpha_smooth * As.transpose() * bs;
            }

            //	Third independent section: POINT CORRESPONDENCE
            //	Construct Ap, bp, A2 and b2
    #pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct Ap and bp
                //------------------------------------------------------------------------
                std::vector<double> bpVector;
                size_t pOffset = 0;

                for (int con_idx = 0; con_idx < constraints.size(); ++con_idx)
                {
                    // size_t vIdx = vertexPair[con_idx].first;
                    // size_t cIdx = vertexPair[con_idx].second;
                    // VertexHandle vh = currentNonRigidMesh.vertex_handle(vIdx);
                    // VertexHandle ch = currentScanMesh.vertex_handle(cIdx);

                    Eigen::Vector3d c,v;
                    v(0) = constraints[con_idx].first.x;
                    v(1) = constraints[con_idx].first.y;
                    v(2) = constraints[con_idx].first.z;
                    c(0) = constraints[con_idx].second.x;
                    c(1) = constraints[con_idx].second.y;
                    c(2) = constraints[con_idx].second.z;

                    double rhs[3] = { 0.0 };

                    for(int con_neigh_idx = 0; 
                        con_neigh_idx < con_node_neigh[con_idx].size();
                        ++ con_neigh_idx)
                    {
                        // actual node index
                        size_t nid = con_node_neigh[con_idx][con_neigh_idx];
                        Eigen::Vector3d n = node_pos[nid];

                        double weight = con_node_weights[con_idx][con_neigh_idx]; // con_node_constraints

                        double vec[3] = { 0.0 };
                        vec[0] = weight * (v[0] - n[0]);
                        vec[1] = weight * (v[1] - n[1]);
                        vec[2] = weight * (v[2] - n[2]);

                        ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nid + 0, vec[0]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nid + 3, vec[1]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nid + 6, vec[2]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 0, 12 * nid + 9, weight));

                        ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nid + 1, vec[0]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nid + 4, vec[1]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nid + 7, vec[2]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 1, 12 * nid + 10, weight));

                        ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nid + 2, vec[0]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nid + 5, vec[1]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nid + 8, vec[2]));
                        ApList.push_back(Eigen::Triplet<double>(pOffset + 2, 12 * nid + 11, weight));

                        rhs[0] -= weight * n[0];
                        rhs[1] -= weight * n[1];
                        rhs[2] -= weight * n[2];
                    }

                    rhs[0] += c[0];
                    rhs[1] += c[1];
                    rhs[2] += c[2];

                    bpVector.push_back(rhs[0]);
                    bpVector.push_back(rhs[1]);
                    bpVector.push_back(rhs[2]);

                    pOffset += 3;
                }

                Ap.resize(pOffset, 12 * node_node_neigh.size());
                Ap.setFromTriplets(ApList.begin(), ApList.end());

                bp.resize(bpVector.size());
                for (int i = 0; i < bpVector.size(); ++i)
                {
                    bp(i) = bpVector[i];
                }

                A2 = params_.m_alpha_point * Eigen::SparseMatrix<double>(Ap.transpose()) * Ap;
                b2 = params_.m_alpha_point * Ap.transpose() * bp;
            }

            //	Fourth independent section: PLANE CORRESPONDENCE
            //	Construct Aq, bq, A3 and b3
            #pragma omp section
            {
                //------------------------------------------------------------------------
                //	Construct Aq and bq
                //------------------------------------------------------------------------
                std::vector<double> bqVector;
                size_t qOffset = 0;

                for (int con_idx = 0; con_idx < constraints.size(); ++con_idx)
                {
                    Eigen::Vector3d c,v, c_n;
                    v(0) = constraints[con_idx].first.x;
                    v(1) = constraints[con_idx].first.y;
                    v(2) = constraints[con_idx].first.z;
                    c(0) = constraints[con_idx].second.x;
                    c(1) = constraints[con_idx].second.y;
                    c(2) = constraints[con_idx].second.z;
                    c_n(0) = constraints[con_idx].first.normal_x;
                    c_n(1) = constraints[con_idx].first.normal_y;
                    c_n(2) = constraints[con_idx].first.normal_z;

                    double rhs = 0.0;

                    for(int con_neigh_idx = 0; 
                        con_neigh_idx < con_node_neigh[con_idx].size();
                        ++ con_neigh_idx)
                    {
                        // actual node index
                        size_t nid = con_node_neigh[con_idx][con_neigh_idx];
                        Eigen::Vector3d n = node_pos[nid];
                        double weight = con_node_weights[con_idx][con_neigh_idx]; // con_node_constraints


                        double vec[3] = { 0.0 };
                        vec[0] = weight * (v[0] - n[0]);
                        vec[1] = weight * (v[1] - n[1]);
                        vec[2] = weight * (v[2] - n[2]);

                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 0, c_n[0] * vec[0]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 1, c_n[1] * vec[0]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 2, c_n[2] * vec[0]));

                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 3, c_n[0] * vec[1]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 4, c_n[1] * vec[1]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 5, c_n[2] * vec[1]));

                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 6, c_n[0] * vec[2]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 7, c_n[1] * vec[2]));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 8, c_n[2] * vec[2]));

                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 9, c_n[0] * weight));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 10, c_n[1] * weight));
                        AqList.push_back(Eigen::Triplet<double>(qOffset, 12 * nid + 11, c_n[2] * weight));

                        rhs -= c_n[0] * weight * n[0];
                        rhs -= c_n[1] * weight * n[1];
                        rhs -= c_n[2] * weight * n[2];
                    }

                    rhs += c_n[0] * c[0];
                    rhs += c_n[1] * c[1];
                    rhs += c_n[2] * c[2];

                    bqVector.push_back(rhs);

                    ++qOffset;
                }

                Aq.resize(qOffset, 12 * node_node_neigh.size());
                Aq.setFromTriplets(AqList.begin(), AqList.end());

                bq.resize(bqVector.size());
                for (int i = 0; i < bqVector.size(); ++i)
                {
                    bq(i) = bqVector[i];
                }

                A3 = params_.m_alpha_plane * Eigen::SparseMatrix<double>(Aq.transpose()) * Aq;
                b3 = params_.m_alpha_plane * Aq.transpose() * bq;
            }
        }

        //------------------------------------------------------------------------
        //	Construct linear equation A * x = b
        //	A = m_alpha_rigid * JrT * Jr + 
        //		m_alpha_smooth * AsT * As +
        //		m_alpha_point * ApT * Ap +
        //		m_alpha_plane * AqT * Aq +
        //
        //	b = - m_alpha_rigid * JrT * fr + 
        //		  m_alpha_smooth * AsT * bs +
        //		  m_alpha_point * ApT * bp +
        //		  m_alpha_plane * AqT * bq +
        //		  m_alpha_rigid * JrT * Jr * affine_vector
        //
        //	x = nextaffine_vector
        //------------------------------------------------------------------------
        Eigen::SparseMatrix<double> A = A0 + A1 + A2 + A3;
        Eigen::VectorXd b = b0 + b1 + b2 + b3;

        // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
        // solver.compute(A);
        // Eigen::VectorXd nextaffine_vector = solver.solve(b);

        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> CholmodSolver;
        CholmodSolver.compute(A);
        Eigen::VectorXd nextaffine_vector = CholmodSolver.solve(b);

        //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> SimplicialLDLTSolver;
        //SimplicialLDLTSolver.compute(A);
        //Eigen::VectorXd nextaffine_vector = SimplicialLDLTSolver.solve(b);

        //------------------------------------------------------------------------
        //	Calculate total energy energyGaussNewton
        //	Break the loop if relative energy change < 1e-6f,
        //	or total iteration is greater than 50
        //------------------------------------------------------------------------
        Eigen::VectorXd RigidVector(12 * node_node_neigh.size());
        Eigen::VectorXd SmoothVector(12 * node_node_neigh.size());
        Eigen::VectorXd PointVector(12 * node_node_neigh.size());
        Eigen::VectorXd PlaneVector(12 * node_node_neigh.size());

        double RigidEnergy = 0.0;
        double SmoothEnergy = 0.0;
        double PointEnergy = 0.0;
        double PlaneEnergy = 0.0;

        //	Parallelize energy calculation
    #pragma omp sections
        {
            //	First section calculates rigid energy
    #pragma omp section
            {
                RigidVector = fr + Jr * (nextaffine_vector - affine_vector.cast<double>());
                RigidEnergy = params_.m_alpha_rigid * RigidVector.transpose() * RigidVector;
            }
            //	Second section calculates smooth energy
    #pragma omp section
            {
                SmoothVector = As * nextaffine_vector - bs;
                SmoothEnergy = params_.m_alpha_smooth * SmoothVector.transpose() * SmoothVector;
            }
            //	Third section calculates point-to-point energy
    #pragma omp section
            {
                PointVector = Ap * nextaffine_vector - bp;
                PointEnergy = params_.m_alpha_point * PointVector.transpose() * PointVector;
            }
            //	Fourth section calculates point-to-plane energy
    #pragma omp section
            {
                PlaneVector = Aq * nextaffine_vector - bq;
                PlaneEnergy = params_.m_alpha_plane * PlaneVector.transpose() * PlaneVector;
            }
        }

        //	Sum up all four energy
        gaussNewtonEnergy = RigidEnergy + SmoothEnergy + PointEnergy + PlaneEnergy;

        // BOOST_LOG_TRIVIAL(debug) 
        // << "\ne_rigid: "<<RigidEnergy
        // << "\ne_smooth: "<<SmoothEnergy
        // << "\ne_point: "<<PointEnergy
        // << "\ne_plane: "<<PlaneEnergy;

        //	Break Gauss-Newton loop if relative energy change is less than 1e-6,
        //	or total iteration is greater than 10
        double gaussNewtonEnergyChange = fabs(gaussNewtonEnergy - gaussNewtonEnergyPrev) / (gaussNewtonEnergyPrev + 1.0);
        // if (verbose)
        // {
            BOOST_LOG_TRIVIAL(debug) 
                << "\tGN iteration " << gaussNewtonIterCnt << ", "
                << "gnEnergy = " << gaussNewtonEnergy 
                << ", gnEnergyChange = " << gaussNewtonEnergyChange;
        // }

        affine_vector = nextaffine_vector;

        if (gaussNewtonEnergyChange < 5e-6 || gaussNewtonIterCnt >= 10) break;
        gaussNewtonEnergyPrev = gaussNewtonEnergy;

        ++gaussNewtonIterCnt;
    }
    BOOST_LOG_TRIVIAL(debug) << "----";

    return gaussNewtonEnergy;

    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////