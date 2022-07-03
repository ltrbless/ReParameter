#include "meshIO.h"
#include "reparam.h"
#include "CLI11.hpp"

int main(int argc, char **argv)
{
	CLI::App app{"Reparam"};

	std::string input_filename;
	int reparam_way = -1;
	bool exportVTK = false;

	app.add_option("-i", input_filename, "input filename. (string, required, supported format: vtk, mesh, pls, obj)")->required()->expected(1, 3);
	app.add_option("-p", reparam_way, "input reparameter wat. (0:tutte)")->required();
	app.add_flag("-k", exportVTK, "Write mesh in VTK format.");

	try
	{
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError &e)
	{
		return app.exit(e);
	}

	size_t input_dotpos = input_filename.find_last_of('.');
	std::string input_postfix = input_filename.substr(input_dotpos + 1, input_filename.length() - input_dotpos - 1);

	Mesh mesh;
	Eigen::MatrixXd V_uv;
	Eigen::MatrixXi T_filled;
	std::string suffix = "";

	if (input_postfix == "vtk")
		MESHIO::readVTK(input_filename, mesh, "surface_id");
	else if(input_postfix == "obj")
		MESHIO::readOBJ(input_filename, mesh);
	else
	{
		std::cout << "file format is not suppport.\n";
		exit(0);
	}

	if(reparam_way != -1)
	{
		switch (reparam_way)
		{
		case 0:
			suffix = ".tutte";
			tiger::buildTuttleParameter(mesh.Vertex, mesh.Topo, V_uv);
			for(int i = 0; i < V_uv.rows(); i++)
			{
				mesh.Vertex(i, 0) = V_uv(i, 0);
				mesh.Vertex(i, 1) = V_uv(i, 1);
				mesh.Vertex(i, 2) = 0;
			}
			break;
		case 1:
			suffix = ".harmonic";
			tiger::buildHarmonicParameter(mesh.Vertex, mesh.Topo, V_uv);
			for(int i = 0; i < V_uv.rows(); i++)
			{
				mesh.Vertex(i, 0) = V_uv(i, 0);
				mesh.Vertex(i, 1) = V_uv(i, 1);
				mesh.Vertex(i, 2) = 0;
			}
			break;
		case 2:
			suffix = ".holeTuttle";
			tiger::buildTuttleHoleParameter(mesh.Vertex, mesh.Topo, V_uv, T_filled);
			mesh.Topo = T_filled;
			mesh.Vertex.resize(V_uv.rows(), 3);
			for(int i = 0; i < V_uv.rows(); i++)
				mesh.Vertex.row(i) << V_uv(i, 0), V_uv(i, 1), 0;
			break;
		default:
			break;
		}
	}

	if (exportVTK)
	{
		std::string output_filename = input_filename.substr(0, input_dotpos) +  suffix +  ".o.vtk";
		MESHIO::writeVTK(output_filename, mesh, "surface_id");
	}





}