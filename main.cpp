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

using namespace std;


class product {
public:
	product(): id_(-1), asin_(""), title_(""), group_(""), avgRating_(-1) {}
	product(int id): id_(id), asin_(""), title_(""), group_(""), avgRating_(-1) {}
	product(int id, const string& asin, const string& title, const string& group, float avg_rating):
		id_(id), asin_(asin), title_(title), group_(group), avgRating_(avg_rating) {}

	int getId() const {return id_;}
	float getRating() const {return avgRating_;}
	const string& getTitle() const {return title_;}
	const string& getAsin() const {return asin_;}
	const string& getGroup() const {return group_;}
	const vector<string>& getSProducts() const {return Simproducts_;}

	void addSProduct(const string& item) {Simproducts_.push_back(item);}
	void addRating(float avg_rating) 	{avgRating_ = avg_rating;}
	void addTitle(const string& title)  {title_ = title;}
	void addAsin(const string& asin)   {asin_ = asin;}
	void addGroup(const string& group)  {group_ = group;}
	void print() const {
		cout << "ID: " << id_ <<  "\tAverage-rating: " << avgRating_ << "\tASIN: " << asin_ << std::endl;
		cout << "Title: " << title_ << "\tGroup: " << group_  << endl;
		cout << "Similar Products: ";
		if (Simproducts_.size() == 0)
			cout << "None";
		for (int i = 0; i < Simproducts_.size(); ++i)
			cout << Simproducts_[i] << " ";
		cout << endl;
	}

private:
	int id_;
	float avgRating_;
	string title_;
	string asin_;
	string group_;
	vector<string> Simproducts_;
};

struct graph {
	int num_verts;
	int num_edges;

	int* degrees;
	int** edges;
	float** edge_weights;
	float* vertex_weights;
};

void read_file(char* filename, int& num_verts, int& num_edges, int*& srcs, int*& dsts,
               std::unordered_map<int, int>& id_Gid, std::map<int, product>& hashToActual) //std::map<int,int>& hashToActual)
{
	ifstream infile;
	string line;
	infile.open(filename);

	num_verts = 0;
	int count = 0, val = 0;
	int alloc_size = 1024;

	srcs = (int*)malloc(alloc_size * sizeof(int));
	dsts = (int*)malloc(alloc_size * sizeof(int));

	while (getline(infile, line)) {
		if (line.c_str()[0] == '%') continue;
		stringstream ss(line);
		if (count + 1 > alloc_size) {
			srcs = (int*)realloc(srcs, alloc_size * 2 * sizeof(int));
			dsts = (int*)realloc(dsts, alloc_size * 2 * sizeof(int));
			alloc_size *= 2;
		}
		ss >> val;
		if (id_Gid.count(val) == 0) {
			id_Gid[val] = num_verts++;
			product tmp(val);
			hashToActual[id_Gid[val]] = tmp;
		}
		srcs[count] = id_Gid[val];

		ss >> val;
		if (id_Gid.count(val) == 0) {
			id_Gid[val] = num_verts++;
			product tmp(val);
			hashToActual[id_Gid[val]] = tmp;
		}
		dsts[count] = id_Gid[val];

		++count;
	}

	num_edges = count;

	infile.close();

	return;
}


void read_Products(char* filename,  std::unordered_map<int, int>& id_Gid, std::map<int, product>& hashToActual,
                   std::unordered_map<string, int>& group_totals) {

	ifstream infile;
	string line, val;
	infile.open(filename);

	std::map<int, product>::iterator itr;

	while (getline(infile, line)) {
		stringstream ss(line);

		ss >> val;
		if (val == "Id:") {
			ss >> val;
			itr = hashToActual.find(id_Gid[atoi(val.c_str())]);
		}
		else if (val == "ASIN:") {
			ss >> val;
			itr->second.addAsin(val);
		}
		else if (val == "title:") {
			ss >> val;
			string title = val;
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
			ss >> val; //throw away # of smiliar
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
		}
	}


	infile.close();
	return;
}

graph* create_graph(int num_verts, int num_edges,
                    int* srcs, int* dsts, float* edge_weights, float* vertex_weights)
{
	graph* g = (graph*)malloc(sizeof(graph));
	g->num_verts = num_verts;
	g->num_edges = num_edges;
	g->degrees = new int[num_verts];
	g->edges = new int*[num_verts];
	g->edge_weights = NULL;
	g->vertex_weights = NULL;
	if (edge_weights != NULL)
		g->edge_weights = new float*[num_verts];
	if (vertex_weights != NULL)
		g->vertex_weights = vertex_weights;

	for (int i = 0; i < num_verts; ++i)
		g->degrees[i] = 0;
	for (int i = 0; i < num_edges; ++i)
		g->degrees[srcs[i]]++;
	for (int i = 0; i < num_edges; ++i)
		g->degrees[dsts[i]]++;
	for (int i = 0; i < num_verts; ++i)
		g->edges[i] = new int[g->degrees[i]];
	if (edge_weights != NULL)
		for (int i = 0; i < num_verts; ++i)
			g->edge_weights[i] = new float[g->degrees[i]];
	for (int i = 0; i < num_verts; ++i)
		g->degrees[i] = 0;
	for (int i = 0; i < num_edges; ++i)
		g->edges[srcs[i]][g->degrees[srcs[i]]++] = dsts[i];
	if (edge_weights != NULL) {
		for (int i = 0; i < num_verts; ++i)
			g->degrees[i] = 0;
		for (int i = 0; i < num_edges; ++i)
			g->edge_weights[srcs[i]][g->degrees[srcs[i]]++] = edge_weights[i];
	}

	return g;
}

void clear_graph(graph*& g)
{
	delete [] g->degrees;
	for (int i = 0; i < g->num_verts; ++i)
		delete [] g->edges[i];
	delete [] g->edges;

	g->num_verts = 0;
	g->num_edges = 0;
}

double eval_density(graph*& g, int* comms)
{
	int* comm_verts = new int[g->num_verts];
	int* int_edges = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		comm_verts[i] = 0;
		int_edges[i] = 0;
	}

	for (int vert = 0; vert < g->num_verts; ++vert) {
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			++comm_verts[comms[vert]];
			if (comms[out] == comms[vert])
				++int_edges[comms[vert]];
		}
	}

	double total_density = 0.0;
	int num_comms = 0;
	for (int i = 0; i < g->num_verts; ++i) {
		if (comm_verts[i] > 1) {
			total_density += (double)int_edges[i] /
			                 ((double)(comm_verts[i] * (comm_verts[i] - 1)) / 2.0);
			++num_comms;
		}
	}

	total_density /= (double)num_comms;

	return total_density;
}


double eval_cutratio(graph*& g, int* comms)
{
	double* num_cuts = new double[g->num_verts];
	double* comm_size = new double[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		num_cuts[i] = 0.0;
		comm_size[i] = 0.0;
	}

	for (int vert = 0; vert < g->num_verts; ++vert) {
		comm_size[comms[vert]] += 1.0;
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			if (comms[out] != comms[vert])
				num_cuts[vert] += 1.0;
		}
	}

	double cutratio = 0.0;
	int num_comms = 0;
	for  (int i = 0; i < g->num_verts; ++i) {
		if (comm_size[i] > 1) {
			++num_comms;
			cutratio += num_cuts[i] / (comm_size[i] * ((double)g->num_verts - comm_size[i]));
		}
	}

	cutratio /= (double)num_comms;

	return cutratio;
}


double eval_conducatance(graph*& g, int* comms)
{
	int* ext_edges = new int[g->num_verts];
	int* int_edges = new int[g->num_verts];
	for (int i = 0; i < g->num_verts; ++i) {
		ext_edges[i] = 0;
		int_edges[i] = 0;
	}

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
			total_conductance += (double)ext_edges[i] / (double)(ext_edges[i] + int_edges[i]);
			++num_comms;
		}
	}

	total_conductance /= (double)num_comms;

	return total_conductance;
}


double eval_modularity(graph*& g, int* comms)
{
	double modularity = 0.0;

	for (int vert = 0; vert < g->num_verts; ++vert) {
		for (int j = 0; j < g->degrees[vert]; ++j) {
			int out = g->edges[vert][j];
			if (comms[out] == comms[vert])
				modularity += (1.0 - (double)(g->degrees[vert] * g->degrees[out]) / (double)(g->num_edges) );
		}
	}

	modularity /= (double)(g->num_edges);

	return modularity;
}

int randomize_queue(int* queue, int size)
{
	for (int i = 0; i < size; ++i) {
		int temp_index = (unsigned)rand() % size;
		int temp_val = queue[temp_index];
		queue[temp_index] = queue[i];
		queue[i] = temp_val;
	}

	return 0;
}

int relabel_comms(int N, int* comms)
{
  int* map = new int[N];
  for (int i = 0; i < N; ++i)
    map[i] = -1;

  int new_num = 0;
  for (int i = 0; i < N; ++i) {
    int comm = comms[i];
    if (map[comm] == -1)
      map[comm] = new_num++;
    comms[i] = map[comm];
  }
  delete [] map;

  return new_num;
}


int* label_prop(graph* g)
{
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
			vector<int> max_labels;
			int max_label = 0;

			for (int j = 0; j < g->degrees[vert]; ++j) {
				int out = g->edges[vert][j];
				int out_label = labels[out];
				if (label_counts.count(out_label) == 0)
					label_counts[out_label] = 1;
				else
					label_counts[out_label] = label_counts[out_label] + 1;

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
			if (max_labels.size() > 1)
				max_index = (unsigned)rand() % max_labels.size();
			if (max_labels.size() > 0 && max_labels[max_index] != labels[vert]) {
				++changes;
				labels[vert] = max_labels[max_index];
			}
		}
	}
	printf("\n");

	return labels;
}



int main(int argc, char** argv)
{
	int* srcs;
	int* dsts;
	int num_verts;
	int num_edges;
	std::unordered_map<int, int> id_Gid;
	std::map<int, product> hashToActual;

	//Book, DVD, Video, Music
	std::unordered_map<string, int> group_totals;

	read_file(argv[1], num_verts, num_edges, srcs, dsts, id_Gid, hashToActual);
	read_Products(argv[2], id_Gid, hashToActual, group_totals);

	graph* g = create_graph(num_verts, num_edges, srcs, dsts, NULL, NULL);
	delete [] srcs;
	delete [] dsts;

	int* comms = label_prop(g);
	printf("Community Structure:\n");
	printf("Density: %lf\n", eval_density(g, comms));
	printf("Cut Ratio: %lf\n", eval_cutratio(g, comms));
	printf("Conductance: %lf\n", eval_conducatance(g, comms));
	printf("Modularity: %lf\n", eval_modularity(g, comms));
	printf("Communities: %d\n", relabel_comms(g->num_verts, comms));


	cout << "There are ";
	for ( const auto& x : group_totals)
		std::cout << x.second << " " << x.first << "s, "; cout << "in the dataset."; cout << endl;

	cout << endl;
	for(std::map<int, product>::iterator itr = hashToActual.begin(); itr != hashToActual.end(); ++itr){
		cout << "Node: " << itr->first << " "; itr->second.print(); cout << endl;
	}

	clear_graph(g);

	return 0;

}
