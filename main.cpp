#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <cfloat>
#include <iostream>
#include <map>
#include "product.h"


struct graph {
	int num_verts;
	int num_edges;

	int* degrees;
	int** edges;
};


//Create source and destination arrays, load important variables
void read_file(char* filename, int& num_verts, int& num_edges, int*& srcs, int*& dsts,
               std::unordered_map<int, int>& IdToGid, std::unordered_map<int, Product>& GidToProd) {

	std::ifstream infile;
	std::string line, sval;
	infile.open(filename);

	num_verts = 0;
	int alloc_size = 1024, count = 0, val = 0;

	srcs = (int*)malloc(alloc_size * sizeof(int));
	dsts = (int*)malloc(alloc_size * sizeof(int));

	while (getline(infile, line)) {
		std::stringstream ss(line);

		//Reallocate arrays to 2x the original size
		if (count + 1 > alloc_size) {
			alloc_size *= 2;
			srcs = (int*)realloc(srcs, alloc_size * sizeof(int));
			dsts = (int*)realloc(dsts, alloc_size * sizeof(int));
		}

		//Load source and destination arrays
		ss >> sval;
		val = atoi(sval.c_str());
		if (val >= 54957)					//Restrict graph size, due to "missing" information
			infile.close();

		// If not in hash table, insert and map true value
		if (IdToGid.count(val) == 0) {
			IdToGid[val] = num_verts++;
			GidToProd[IdToGid[val]] = Product(val);
		}
		srcs[count] = IdToGid[val];

		ss >> sval;
		val = atoi(sval.c_str());
		if (IdToGid.count(val) == 0) {
			IdToGid[val] = num_verts++;
			GidToProd[IdToGid[val]] = Product(val);
		}
		dsts[count] = IdToGid[val];

		++count;
	}

	num_edges = count;

	infile.close();

	return;
}

//Add product information for each node
void read_Products(char* filename,  std::unordered_map<int, int>& IdToGid, std::unordered_map<int, Product>& GidToProd,
                   std::unordered_map<std::string, int>& group_totals) {

	std::ifstream infile;
	std::string line, val;
	infile.open(filename);

	std::unordered_map<int, Product>::iterator itr;

	while (getline(infile, line)) {
		std::stringstream ss(line);

		ss >> val;
		if (val == "Id:") {
			ss >> val;
			itr = GidToProd.find(IdToGid[atoi(val.c_str())]);
		}
		else if (val == "ASIN:") {
			ss >> val;
			itr->second.addAsin(val);
		}
		else if (val == "title:") {
			ss >> val;
			std::string title = val;
			while (ss >> val) {
				title += " " + val;
			}
			itr->second.addTitle(title);
		}
		else if (val == "group:") {
			ss >> val;
			itr->second.addGroup(val);
			group_totals[val]++;
		}
		else if (val == "similar:") {
			ss >> val; 									//throw away # of smiliar
			while (ss >> val)
				itr->second.addSProduct(val);
		}
		else if (val == "rating:") {
			ss >> val;
			itr->second.addRating(atof(val.c_str()));
		}
		else if (val == "discontinued") {
			itr->second.addTitle(val + " product");
			itr->second.addGroup(val + " product");
			group_totals[val]++;
		}
	}


	infile.close();
	return;
}

graph* create_graph(int num_verts, int num_edges, int* srcs, int* dsts) {
	graph* g = (graph*)malloc(sizeof(graph));
	g->num_verts = num_verts;
	g->num_edges = num_edges;
	g->degrees = new int[num_verts];
	g->edges = new int*[num_verts];


	for (int i = 0; i < num_verts; ++i)
		g->degrees[i] = 0;

	//Find total degrees for all vertices
	for (int i = 0; i < num_edges; ++i)
		g->degrees[srcs[i]]++;

	// Create adjacency array from degrees
	for (int i = 0; i < num_verts; ++i)
		g->edges[i] = new int[g->degrees[i]];
	for (int i = 0; i < num_verts; ++i)
		g->degrees[i] = 0;

	//Fill out adjacency array
	for (int i = 0; i < num_edges; ++i)
		g->edges[srcs[i]][g->degrees[srcs[i]]++] = dsts[i];

	return g;
}


void clear_graph(graph*& g) {
	delete [] g->degrees;
	for (int i = 0; i < g->num_verts; ++i)
		delete [] g->edges[i];
	delete [] g->edges;

	g->num_verts = 0;
	g->num_edges = 0;
}


double eval_density(graph*& g, int* comms, int& num_comms) {
	int* comm_verts = new int[g->num_verts];
	int* int_edges = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		comm_verts[i] = 0;
		int_edges[i] = 0;
	}

	//count occurrence of labels && check if neighbors have the same label
	for (int vert = 0; vert < g->num_verts; ++vert) {
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			++comm_verts[comms[vert]];
			if (comms[out] == comms[vert])
				++int_edges[comms[vert]];
		}
	}

	//Calculate density for each community
	double total_density = 0.0;
	num_comms = 0;
	for (int i = 0; i < g->num_verts; ++i) {
		if (comm_verts[i] > 1) {
			total_density += int_edges[i] /
			                 ( (double)(comm_verts[i] * (comm_verts[i] - 1)) / 2.0);
			++num_comms;
		}
	}

	total_density /= num_comms;

	delete [] comm_verts;
	delete [] int_edges;

	return total_density;
}


double eval_cutratio(graph*& g, int* comms) {
	int* num_cuts = new int[g->num_verts];
	int* comm_size = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		num_cuts[i] = 0;
		comm_size[i] = 0;
	}

	// Count number of edge cuts we have to preform
	for (int vert = 0; vert < g->num_verts; ++vert) {
		comm_size[comms[vert]] += 1;
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			if (comms[out] != comms[vert])
				num_cuts[vert] += 1;
		}
	}

	//Calculate cutratio for each community
	double cutratio = 0.0;
	int num_comms = 0;
	for  (int i = 0; i < g->num_verts; ++i) {
		if (comm_size[i] > 1) {
			++num_comms;
			cutratio += num_cuts[i] / (double)(comm_size[i] * (g->num_verts - comm_size[i]));
		}
	}

	cutratio /= num_comms;

	delete [] num_cuts;
	delete [] comm_size;

	return cutratio;
}


double eval_conducatance(graph*& g, int* comms) {
	int* ext_edges = new int[g->num_verts];
	int* int_edges = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		ext_edges[i] = 0;
		int_edges[i] = 0;
	}

	//Check if our neighbors are within our commmunity
	for (int vert = 0; vert < g->num_verts; ++vert) {
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			if (comms[out] != comms[vert])
				++ext_edges[comms[vert]];
			else
				++int_edges[comms[vert]];
		}
	}

	double total_conductance = 0.0;
	int num_comms = 0;
	for (int i = 0; i < g->num_verts; ++i) {
		if (ext_edges[i] > 0 || int_edges[i] > 0) {
			total_conductance += ext_edges[i] / (double)(ext_edges[i] + int_edges[i]);
			++num_comms;
		}
	}

	total_conductance /= num_comms;

	delete [] ext_edges;
	delete [] int_edges;

	return total_conductance;
}


double eval_modularity(graph*& g, int* comms) {
	double modularity = 0.0;

	//Check if our neighbors are within our commmunity
	for (int vert = 0; vert < g->num_verts; ++vert) {
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			if (comms[out] == comms[vert])
				modularity += (1.0 - (double)(g->degrees[vert] * g->degrees[out]) / (g->num_edges) );
		}
	}

	modularity /= (g->num_edges);

	return modularity;
}


int randomize_queue(int* queue, int size) {
	for (int i = 0; i < size; ++i) {
		srand(time(NULL));
		int temp_index = (unsigned)rand() % size;
		int temp_val = queue[temp_index];
		queue[temp_index] = queue[i];
		queue[i] = temp_val;
	}

	return 0;
}


int* label_prop(graph* g) {
	int max_iter = 20;

	int* labels = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i)
		labels[i] = i;

	int* queue = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i)
		queue[i] = i;

	int changes = 1;
	int num_iter = 0;
	while (changes && num_iter < max_iter)
	{
		changes = 0;
		randomize_queue(queue, g->num_verts);
		++num_iter;
		for (int j = 0; j < g->num_verts; ++j) {
			int vert = queue[j];

			std::unordered_map<int, int> label_counts;
			std::vector<int> max_labels;
			int max_label = 0;

			//For all of vert's neighbors
			for (int j = 0; j < g->degrees[vert]; ++j) {
				int out = g->edges[vert][j];
				int out_label = labels[out];

				//Check if we've seen the label before
				if (label_counts.count(out_label) == 0)
					label_counts[out_label] = 1;
				else
					label_counts[out_label] = label_counts[out_label] + 1;

				//Check if this is the most seen label
				if (label_counts[out_label] > max_label) {
					max_labels.clear();
					max_labels.push_back(out_label);
					max_label = label_counts[out_label];
				}
				else if (label_counts[out_label] == max_label) {
					max_labels.push_back(out_label);
				}
			}

			int max_index = 0;

			//If mutliple -> any label can replace
			if (max_labels.size() > 1)
				max_index = (unsigned)rand() % max_labels.size();

			//If single, highest label will overwrite. if mulitple then all labels have 1/max_label.size() prob of overwriting
			//Will converage to highest label over multiple iterations
			if (max_labels.size() > 0 && max_labels[max_index] != labels[vert]) {
				++changes;
				labels[vert] = max_labels[max_index];
			}
		}
	}

	delete [] queue;
	return labels;
}


void sample_community(graph* g, int* comms, const std::unordered_map<int, Product>& GidToProd) {
	srand(time(NULL));
	int ran_index = (unsigned)rand() % g->num_verts;
	int label = comms[ran_index];

	//Collect all products within the sampled community
	std::vector<Product> community_vertexs;
	for (int vert = 0; vert < g->num_verts; ++vert)
		if (comms[vert] == label)
			community_vertexs.push_back(GidToProd.at(vert));

	std::unordered_map<std::string, int> groups;
	float rating = 0.0;

	int missingProducts = 0;
	//Printing sampled community
	for (int i = 0 ; i < community_vertexs.size(); ++i) {
		community_vertexs[i].print();
		groups[community_vertexs[i].getGroup()]++;
		//To account for missing data - will get true rating
		if(community_vertexs[i].getRating() == -1){
			missingProducts++;
			continue;
		}
		rating += community_vertexs[i].getRating();

	}

	//Find largest group within the sampled commmunity
	std::string largest_group;
	int group_size = 0;
	for ( const auto& itr : groups) {
		if (itr.second > group_size) {
			group_size = itr.second;
			largest_group = itr.first;
		}
	}

	if (largest_group.empty()) largest_group = "Unknown item";
	printf("\nSampled community size: %lu\n", community_vertexs.size());
	printf("Community is primarily made up of: %s's, with %d of them.\n", largest_group.c_str(), group_size);
	printf("The average rating for this community is %.2f\n", rating / (double)(community_vertexs.size() - missingProducts) );

}

void print_stats(graph* g, int* comms, const std::unordered_map<std::string, int>& group_totals, int num_verts) {
	int total_comms = 0;

	printf("\nTotal Community Structure:\n");
	printf("Density: %lf\n", eval_density(g, comms, total_comms));
	printf("Cut Ratio: %lf\n", eval_cutratio(g, comms));
	printf("Conductance: %lf\n", eval_conducatance(g, comms));
	printf("Modularity: %lf\n", eval_modularity(g, comms)); \
	printf("Total Communities: %d\nAverage Community size: %d\n", total_comms, num_verts / total_comms);

	std::cout << "\nThere are ";
	for ( const auto& x : group_totals)
		printf("%d %s's, ", x.second, x.first.c_str());

	std::cout << "in the dataset." << std::endl;
	std::cout << "Missing information for 493,595 items." << std::endl;
}


//  Prints out mapping with each nodes metadata
void print_vertex(const std::unordered_map<int, Product>& GidToProd) {
	std::cout << std::endl;
	for (const auto& itr : GidToProd) {
		std::cout << "Node: " << itr.first << " "; itr.second.print(); std::cout << std::endl;
	}
}


int main(int argc, char** argv) {
	int* srcs;
	int* dsts;
	int num_verts;
	int num_edges;

	//Mapings
	std::unordered_map<int, int> IdToGid;
	std::unordered_map<int, Product> GidToProd;

	//Book, DVD, Video, Music
	std::unordered_map<std::string, int> group_totals;

	read_file(argv[1], num_verts, num_edges, srcs, dsts, IdToGid, GidToProd);
	read_Products(argv[2], IdToGid, GidToProd, group_totals);
	graph* g = create_graph(num_verts, num_edges, srcs, dsts);

	delete [] srcs;
	delete [] dsts;

	int* comms = label_prop(g);

	print_stats(g, comms, group_totals, num_verts);

	std::cout << "\nSampled Community:" << std::endl;
	sample_community(g, comms, GidToProd);

	//Option to print all vertices
	char ans;
	std::cout << "\nWould you like to see the data regarding each vertex?(y/n): ";
	std::cin >> ans;
	if (ans == 'y')
		print_vertex(GidToProd);

	clear_graph(g);

	return 0;

}







