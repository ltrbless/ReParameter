#include "getloop.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <fstream>
#include <set>


namespace tiger
{

void dfs_get_loop2(int cur, int pre, std::vector<bool>& vis, std::vector<std::vector<int>>& G, std::vector<int>& path, std::vector<std::vector<int>>& loop_lst)
{
	// std::cout << cur << " " << vis[cur] << "\n";
	if(vis[cur])
	{
		std::vector<int> tmp;
		for(int i = path.size() - 2; i >= 0; i--)
		{
			// std::cout << path[i] << " * " << cur << "\n";
			if(path[i] != cur)
			{
				tmp.push_back(path[i]);
			}
			else
			{
				tmp.push_back(path[i]);
				break;
			}
		}
		loop_lst.push_back(tmp);
		return ;
	}

	vis[cur] = 1;
	for(int i = 0; i < G[cur].size(); i++)
	{
		if(G[cur][i] == pre) continue;
		path.push_back(G[cur][i]);
		dfs_get_loop2(G[cur][i], cur, vis, G, path, loop_lst);
		path.pop_back();
	}
	vis[cur] = 0;

}

void boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& bnd)
{
	// 检测是否有重复点
	double eps = 1e-3;
	Eigen::MatrixXd V_tmp;
	Eigen::MatrixXi F_tmp;
	Eigen::VectorXi ind;

	V_tmp = V;
	F_tmp = F;

	bool b_dulplicated_point = false;
	if(b_dulplicated_point)
	{
		for(int i = 0; i < V.rows(); i++)
		{
			for(int j = i + 1; j < V.rows(); j++)
			if( (V.row(i) - V.row(j)).norm() < eps )
			{
				std::cout << "Warning : boundary loop have dulplicate point.\n";
				break;
			}
		}
	}

	// 找到所有的边界边，并建图
	std::vector<std::vector<int>> G(V.rows(), std::vector<int>());
	std::map<std::array<int, 2> , int> mp;
	std::vector<bool> vis(V.rows(), 0);
	std::vector<std::vector<int>> loop_lst;
	std::vector<int> path;
	std::set<int> num_p_set;


	for(int i = 0; i < F.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3 );
			std::array<int, 2> tmp = {std::min(a, b), std::max(a, b)};
			mp[tmp]++;
		}
	}

	for(auto iter = mp.begin(); iter != mp.end(); iter++)
	{
		if( iter->second == 1 )
		{
			G[iter->first[0]].push_back(iter->first[1]);
			G[iter->first[1]].push_back(iter->first[0]);
			// std::cout << iter->first[0] << " " << iter->first[1] << '\n';
			num_p_set.insert(iter->first[0]);
			num_p_set.insert(iter->first[1]);
		}
	}

	int start = *num_p_set.begin();
	path.push_back(start);
	dfs_get_loop2(start, -1, vis, G, path, loop_lst);

	for(int i = 0; i < loop_lst.size(); i++)
	{
		std::cout << "#BK The " << i << "th loop number is : " << loop_lst[i].size() << '\n';
	}

	sort(loop_lst.begin(), loop_lst.end(), [&V](std::vector<int>& a, std::vector<int>& b){
		double len_a = 0.0;
		double len_b = 0.0;

		for(int i = 0; i < a.size(); i++)
		{
			len_a += (V.row(i) - V.row((i + 1) % a.size())).norm();
		}

		for(int i = 0; i < b.size(); i++)
		{
			len_b += (V.row(i) - V.row((i + 1) % b.size())).norm();
		}
		
		return len_a > len_b;
	});	


	int s = loop_lst[0][0];
	int e = loop_lst[0][1];
	int ok = 0;
	for(int i = 0; i < F.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3 );
			if(a == s && b == e)
			{
				ok = 1;
			}
		}
		if(ok) break;
	}

	if(!ok)
	{
		std::reverse(loop_lst[0].begin(), loop_lst[0].end());
	}

	bnd.resize(loop_lst[0].size(), 1);
	for(int i = 0; i < loop_lst[0].size(); i++)
	{
		bnd(i) = loop_lst[0][i];
	}

}

void all_boundary_loop_by_dfs2(Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Eigen::VectorXi>& bnd_vec)
{
	// 检测是否有重复点
	double eps = 1e-3;
	Eigen::MatrixXd V_tmp;
	Eigen::MatrixXi F_tmp;
	Eigen::VectorXi ind;

	V_tmp = V;
	F_tmp = F;

	bool b_dulplicated_point = false;
	if(b_dulplicated_point)
	{
		for(int i = 0; i < V.rows(); i++)
		{
			for(int j = i + 1; j < V.rows(); j++)
			if( (V.row(i) - V.row(j)).norm() < eps )
			{
				std::cout << "Warning : boundary loop have dulplicate point.\n";
				break;
			}
		}
	}

	// 找到所有的边界边，并建图
	std::vector<std::vector<int>> G(V.rows(), std::vector<int>());
	std::map<std::array<int, 2> , int> mp;
	std::vector<bool> vis(V.rows(), 0);
	std::vector<std::vector<int>> loop_lst;
	std::vector<int> path;
	std::set<int> num_p_set;


	for(int i = 0; i < F.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3 );
			std::array<int, 2> tmp = {std::min(a, b), std::max(a, b)};
			mp[tmp]++;
		}
	}

	for(auto iter = mp.begin(); iter != mp.end(); iter++)
	{
		if( iter->second == 1 )
		{
			G[iter->first[0]].push_back(iter->first[1]);
			G[iter->first[1]].push_back(iter->first[0]);
			// std::cout << iter->first[0] << " " << iter->first[1] << '\n';
			num_p_set.insert(iter->first[0]);
			num_p_set.insert(iter->first[1]);
		}
	}

	int start = *num_p_set.begin();
	path.push_back(start);
	dfs_get_loop2(start, -1, vis, G, path, loop_lst);

	for(int i = 0; i < loop_lst.size(); i++)
	{
		std::cout << "#BK The " << i << "th loop number is : " << loop_lst[i].size() << '\n';
	}

	sort(loop_lst.begin(), loop_lst.end(), [&V](std::vector<int>& a, std::vector<int>& b){
		double len_a = 0.0;
		double len_b = 0.0;

		for(int i = 0; i < a.size(); i++)
		{
			len_a += (V.row(i) - V.row((i + 1) % a.size())).norm();
		}

		for(int i = 0; i < b.size(); i++)
		{
			len_b += (V.row(i) - V.row((i + 1) % b.size())).norm();
		}
		
		return len_a > len_b;
	});	

  // new loop
  std::unordered_map<int, bool> mploop;
  std::vector<std::vector<int>> new_loop;
  
  for(int i = 0; i < loop_lst.size(); i++)
  {
    int num = 0;
    for(int j = 0; j < loop_lst[i].size(); ++j)
    {
      if(mploop[loop_lst[i][j]]) num++;
    }
    if(num < 3)
    {
      new_loop.push_back(loop_lst[i]);
      for(int j = 0; j < loop_lst[i].size(); ++j) mploop[loop_lst[i][j]] = 1;
    }
  }
  loop_lst.clear();
  loop_lst = new_loop;

	int s = loop_lst[0][0];
	int e = loop_lst[0][1];
	int ok = 0;
	for(int i = 0; i < F.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			int a = F(i, j);
			int b = F(i, (j + 1) % 3 );
			if(a == s && b == e)
			{
				ok = 1;
			}
		}
		if(ok) break;
	}

	if(!ok)
	{
		std::reverse(loop_lst[0].begin(), loop_lst[0].end());
	}

	// for(int i = 0; i < loop_lst.size(); ++i)
	// {
	// 	for(int j = 0; j < loop_lst[i].size(); ++j)
	// 	{
	// 		std::cout << loop_lst[i][j] << "  ";
	// 	}
	// 	std::cout << "\n";
	// }

  for(int i = 0; i < loop_lst.size(); ++i)
  {
    Eigen::MatrixXi bnd;
    bnd.resize(loop_lst[i].size(), 1);
    for(int j = 0; j < loop_lst[i].size(); ++j)
    {
      bnd(j) = loop_lst[i][j];
    }
    bnd_vec.push_back(bnd);
  }

}

};