
#include "reparam.h"
#include "getloop.h"
#include "meshIO.h"
#include "kdtree_m.h"

#include <igl/bfs_orient.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/harmonic.h>

#include <set>
#include <unordered_map>
#include <cmath>

namespace tiger
{

void buildTuttleHoleParameter( Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv, Eigen::MatrixXi& T_filled )
{
	Eigen::MatrixXi C;
  std::vector<Eigen::VectorXi> bnd_vec;

	igl::bfs_orient(T_3d, T_3d, C);


  // remove duplicate point
  double minn_size = std::numeric_limits<double>::max();

  for(int i = 0; i < T_3d.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      int a = T_3d(i, j);
      int b = T_3d(i, (j + 1) % 3);
      minn_size = std::min(minn_size, (V_3d.row(a) - V_3d.row(b)).norm() );
    }
  }
  minn_size = sqrt(minn_size);
  Eigen::MatrixXd tmpV = V_3d;
  Eigen::MatrixXi tmpT = T_3d;
  Eigen::MatrixXi SVI, SVJ;
  igl::remove_duplicate_vertices(tmpV, tmpT, minn_size / 1000.0, V_3d, SVI, SVJ, T_3d);

  // removeDulplicatePoint(V_3d, T_3d, minn_size / 1000.0);

  all_boundary_loop_by_dfs2(V_3d, T_3d, bnd_vec); // mine
    
  // fill hole  
  Eigen::MatrixXi& F_filled = T_filled;
  int n_filled_faces = 0;
  int num_holes = bnd_vec.size() - 1;
  int real_F_num = T_3d.rows();
  const int V_rows = T_3d.maxCoeff()+1;

  for (int i = 0; i < num_holes; i++)
    n_filled_faces += bnd_vec[i + 1].size();
  F_filled.resize(n_filled_faces + real_F_num, 3);
  F_filled.topRows(real_F_num) = T_3d;

  int new_vert_id = V_rows;
  int new_face_id = real_F_num;

  for (int i = 0; i < num_holes; i++, new_vert_id++)
  {
    int cur_bnd_size = bnd_vec[i + 1].size();
    int it = 0;
    int back = bnd_vec[i + 1].size() - 1;
    F_filled.row(new_face_id++) << bnd_vec[i + 1][it], bnd_vec[i + 1][back], new_vert_id;
    // std::cout << bnd_vec[i + 1][it] << "  " <<  bnd_vec[i + 1][back] << "  " <<  new_vert_id << "\n";
    while (it != back)
    {
      F_filled.row(new_face_id++)
          << bnd_vec[i + 1][(it + 1)],
          bnd_vec[i + 1][(it)], new_vert_id;
    // std::cout << bnd_vec[i + 1][it] << "  " <<  bnd_vec[i + 1][it + 1] << "  " <<  new_vert_id << "\n";
      it++;
    }
  }
  assert(new_face_id == F_filled.rows());
  assert(new_vert_id == V_rows + num_holes);

  int V_N = new_vert_id;
  // int V_N = V_3d.rows();


  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V_3d, bnd_vec[0], bnd_uv);


  // build adj relation
  std::vector<std::set<int>> adj(V_N);
  for(int i = 0; i < F_filled.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      int a = F_filled(i, j);
      int b = F_filled(i, (j + 1) % 3);
      // std::cout << a << "  " << b << "\n";
      adj[a].insert(b);
      adj[b].insert(a);
    }
  }

  typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SMatrix;

	std::vector<T> tripletlist;
	Eigen::VectorXd bu = Eigen::VectorXd::Zero(V_N);
	Eigen::VectorXd bv = Eigen::VectorXd::Zero(V_N);

  std::unordered_map<int, bool> is_boundary;
  std::unordered_map<int, Eigen::Vector2d> mp_bnd;
  for(int i = 0; i < bnd_vec[0].rows(); i++)
  {
    is_boundary[bnd_vec[0](i, 0)] = 1;
    mp_bnd[bnd_vec[0](i, 0)] = bnd_uv.row(i);
  }

  for(int i = 0; i < V_N; ++i)
  {
    if(is_boundary[i])
    {
      tripletlist.emplace_back(i, i, 1);
      bu(i) = mp_bnd[i].x();
      bv(i) = mp_bnd[i].y();
    }
    else{
      for(auto& t : adj[i])
      {
        tripletlist.emplace_back(i, t, -1);
      }
      tripletlist.emplace_back(i, i, adj[i].size());
    }
  }

  SMatrix coff(V_N, V_N);
  coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
  Eigen::SparseLU<SMatrix> solver;
  solver.compute(coff);

  Eigen::VectorXd xu = solver.solve(bu);
  Eigen::VectorXd xv = solver.solve(bv);

  V_uv.resize(V_N, 2);

  V_uv.col(0) = xu;
  V_uv.col(1) = xv;

  std::cout << "Build fill hole tuttle parameter is success.\n";
}
  

void buildTuttleParameter( Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv )
{
	Eigen::MatrixXi C;
  Eigen::VectorXi bnd;

	igl::bfs_orient(T_3d, T_3d, C);


  // remove duplicate point
  double minn_size = std::numeric_limits<double>::max();

  for(int i = 0; i < T_3d.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      int a = T_3d(i, j);
      int b = T_3d(i, (j + 1) % 3);
      minn_size = std::min(minn_size, (V_3d.row(a) - V_3d.row(b)).norm() );
    }
  }
  minn_size = sqrt(minn_size);
  Eigen::MatrixXd tmpV = V_3d;
  Eigen::MatrixXi tmpT = T_3d;
  Eigen::MatrixXi SVI, SVJ;
  igl::remove_duplicate_vertices(tmpV, tmpT, minn_size / 1000.0, V_3d, SVI, SVJ, T_3d);

  // removeDulplicatePoint(V_3d, T_3d, minn_size / 1000.0);
  int V_N = V_3d.rows();

  boundary_loop_by_dfs2(V_3d, T_3d, bnd); // mine
    
  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V_3d, bnd, bnd_uv);


  // build adj relation
  std::vector<std::vector<int>> adj(V_N);
  for(int i = 0; i < T_3d.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      int a = T_3d(i, j);
      int b = T_3d(i, (j + 1) % 3);
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
  }

  typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SMatrix;

	std::vector<T> tripletlist;
	Eigen::VectorXd bu = Eigen::VectorXd::Zero(V_N);
	Eigen::VectorXd bv = Eigen::VectorXd::Zero(V_N);

  std::unordered_map<int, bool> is_boundary;
  std::unordered_map<int, Eigen::Vector2d> mp_bnd;
  for(int i = 0; i < bnd.rows(); i++)
  {
    is_boundary[bnd(i, 0)] = 1;
    mp_bnd[bnd(i, 0)] = bnd_uv.row(i);
  }

  for(int i = 0; i < V_N; ++i)
  {
    if(is_boundary[i])
    {
      tripletlist.emplace_back(i, i, 1);
      bu(i) = mp_bnd[i].x();
      bv(i) = mp_bnd[i].y();
    }
    else{
      for(int j = 0; j < adj[i].size(); j++)
      {
        tripletlist.emplace_back(i, adj[i][j], -1);
      }
      tripletlist.emplace_back(i, i, adj[i].size());
    }
  }

  SMatrix coff(V_N, V_N);
  coff.setFromTriplets(tripletlist.begin(), tripletlist.end());
  Eigen::SparseLU<SMatrix> solver;
  solver.compute(coff);

  Eigen::VectorXd xu = solver.solve(bu);
  Eigen::VectorXd xv = solver.solve(bv);

  V_uv.resize(V_N, 2);

  V_uv.col(0) = xu;
  V_uv.col(1) = xv;

  std::cout << "Build tuttle parameter is success.\n";
}


void buildHarmonicParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv)
{
  Eigen::VectorXi bnd;

  // remove duplicate point
  double minn_size = std::numeric_limits<double>::max();

  for(int i = 0; i < T_3d.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      int a = T_3d(i, j);
      int b = T_3d(i, (j + 1) % 3);
      minn_size = std::min(minn_size, (V_3d.row(a) - V_3d.row(b)).norm() );
    }
  }
  minn_size = sqrt(minn_size);
  Eigen::MatrixXd tmpV = V_3d;
  Eigen::MatrixXi tmpT = T_3d;
  Eigen::MatrixXi SVI, SVJ;
  igl::remove_duplicate_vertices(tmpV, tmpT, minn_size / 1000.0, V_3d, SVI, SVJ, T_3d);

  boundary_loop_by_dfs2(V_3d, T_3d, bnd); // mine

  // Map the boundary to a circle, preserving edge proportions
  Eigen::MatrixXd bnd_uv;
  igl::map_vertices_to_circle(V_3d, bnd, bnd_uv);

  igl::harmonic(V_3d, T_3d, bnd, bnd_uv, 1, V_uv);

  // Scale UV to make the texture more clear
  // V_uv *= 10000;

  std::cout << "Build harmonic parameter is success.\n";
}

void removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& T, double eps)
{
    std::vector<vec3d> V_out_vec;
    std::vector<int> J(V.rows());
    KdTreeM kdtree(3);
    int baseid = 0;
    for(int i = 0; i < V.rows(); i++)
    {
        vec3d curp = vec3d( V(i, 0), V(i, 1), V(i, 2) );
        std::vector<int> ansvec = kdtree.Query3DNodeByDistance(curp, eps);
        if(ansvec.size() > 0)
        {
            J[i] = ansvec[0];
            continue;
        }
        kdtree.Insert3DNode(curp, baseid);
        V_out_vec.push_back(curp);
        J[i] = baseid;
        baseid++;
    }

    V.resize(V_out_vec.size(), 3);
    for(int i = 0; i < V_out_vec.size(); i++)
    {
        V.row(i) = V_out_vec[i];
    }

    // std::cout << T.rows() << '\n';
    std::vector<std::array<int, 3>> F_out;
    for(int i = 0; i < T.rows(); i++) {
        std::array<int, 3> facet_out;
        for(int j = 0; j < 3; j++) {
            facet_out[j] = J[T(i, j)]; 
        }
        if(facet_out[0] == facet_out[1] || facet_out[1] == facet_out[2] || facet_out[0] == facet_out[2]) {
            continue;
        }
        F_out.push_back(facet_out);
    }
    // std::cout << F_out.size() << '\n';
    T.resize(F_out.size(), 3);
    for(int i = 0; i < F_out.size(); i++) {
        for(int j = 0; j < 3; j++) {
            T(i, j) = F_out[i][j];
        }
    }
}

};