#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <stack>

using namespace std;

#define RANDOM_LIMIT 9
#define n_cities 3
#define hometown 0 
#define NC -1

//en la representación de un digrafo cada fila representará los costos de ida y cada columna los costos de vuelta

class Matrix{
	public:
		int n;
		int m;
		vector <vector <double>> matrix;

		Matrix(int n1, int m1){ //mi constructor llena la matriz con valores aleatorios simulando distancia entre ciudades de un grafo
			n = n1; 
			m = m1;
			vector <vector <double>> matrix1(n);
			matrix = matrix1;
			for(int i = 0; i < n; i++){
				vector <double> row(m);
				for(int j = 0; j < m; j++){
					if(i != j) row[j] = (double)(1+rand()%(RANDOM_LIMIT));
					else row[j] = 0.0;
				}
				matrix[i] = row;
			}
		}

		void print(){
			for(int i = 0; i < n; i++){
				for(int j = 0; j < m; j++){
					cout<<matrix[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
		}
};

class Node{
	public:
		int data;
		vector <int> history;
		double cost;
		int counter; //hasta no alcanzar el máximo número de ciudades se construiran recursivamente más niveles
		vector <Node*> children;
		int before;
		Node *bfr; 
		Node(){
			cost = 0.0;
			before = -1;
		}
		void buildTree(Matrix* &digraph, int hometown1, int counter1){
			//data = hometown1;
			cout<<"hometown1 "<<hometown1<<endl;
			vector <double> ways = digraph->matrix[hometown1];
			Node *nc = new Node();
			nc->data = data;
			nc->history = history;
			nc->bfr = bfr;
			nc->cost = cost;
			nc->counter = counter;
			nc->children = children;
			if(counter1 == 0){
				data = hometown1;
				history.push_back(hometown1);
				nc->data = data;
				cout<<"homet "<<nc->data<<endl;
				nc->history = history;
				for(int i = 0; i < n_cities; i++){
					if(ways[i] != 0){
						Node *n1 = new Node();
						n1->data = i;
						n1->cost += ways[i];
						//n1->history.push_back(hometown1);
						n1->before = hometown1;
						children.push_back(n1);
						nc->children = children;
						n1->bfr = nc;
						counter = counter1;
						
					}
				}
			}else{
				for(int i = 0; i < bfr->history.size(); i++){
					history.push_back(bfr->history[i]);
				}
				//history = bfr->history;
				history.push_back(hometown1);
				nc->history = history;
				vector <int> tmp;
				for(int i = 0; i < n_cities; i++){
					tmp.push_back(0);
				}
				for(int i = 0; i < history.size(); i++){
					tmp[history[i]]+=1;
				}
				cout<<"ht "<<hometown1<<endl;
				for(int i = 0; i < n_cities; i++){
					cout<<"ss "<<tmp[i]<<endl;
				}
				for(int i = 0; i < n_cities; i++){
					if(tmp[i] == 0){
						cout<<"entro"<<endl;
						Node *n1 = new Node();
						n1->data = i;
						cout<<"aca"<<n1->data<<endl;
						n1->cost += ways[i];
						n1->bfr = nc;
						//n1->history = history;
						//n1->history.push_back(i);
						n1->before = hometown1;
						counter = counter1;
						children.push_back(n1);
					}
				}
				

			}
				/*
				for(int i = 0; i < n_cities; i++){
					for(int j = 0; j < history.size(); j++){
						if(history[j] != i){ 
							cout<<"entro"<<endl;
							Node *n1 = new Node();
							n1->data = i;
							n1->cost += ways[i];
							for(int k = 0; k < history.size(); k++){
								n1->history.push_back(history[k]);
							}
							n1->history.push_back(i);
							n1->before = hometown1;
							counter = counter1;
							children.push_back(n1);
						}
					}	
				}
				cout<<endl;
			}
			*/
			
			counter1 += 1;
			cout<<"numero de hijos "<<children.size()<<endl;
			if(children.size() > 0){
				for(int i = 0; i < children.size(); i++){
					children[i]->buildTree(digraph, children[i]->data, counter1);	
				}
			}
		}
		void print(){
			cout<<"("<<data<<") con padre "<<before<<endl;
			if(!children.empty()){
				for(int i = 0; i < children.size(); i++){
					children[i]->print();
				}
			}
			
		}
		//verifica si ya ha sido visitada una ciudad y si conviene tomar una ruta con un mejor costo
		bool feasible(vector <int> tour, int city, int totalCost){
			if(find(tour.begin(), tour.end(), city) != tour.end()){
				return false;
			}else{
				int newTotalCost = totalCost + cost;
				if(newTotalCost > totalCost)
					return false;
				else 
					return true;
			}
		}
		//recursiva
		void depth_first_search(vector <int> &bestTour, vector <int> &tour, int totalCost){
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				for(int i = 0; i < children.size(); i++){
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
		}
		//no recursiva
		void depth_first_search2(stack <int> S, vector <int> &tour){
			for(int i = n_cities-1; i >= 0; i--){
				S.push(i);
			}
			int totalCost = 0;
			while(!S.empty()){
				int city = S.top();
				if(city == NC){
					tour.pop_back();
				}
				else{
					tour.push_back(data);
					if(tour.size() == n_cities){
						if(feasible(tour, city, totalCost)){
							tour.push_back(data);
						}
						tour.pop_back();

					}else{
						S.push(NC);
						for(int i = n_cities-1; i >= 0; i++){
							if(feasible(tour, city, totalCost))
								S.push(city);
						}
					}
				}
			}
		}

		//paralelizada 1
		void depth_first_searchP1(vector <int> &bestTour, vector <int> &tour, int totalCost){
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				#pragma omp parallel for
				for(int i = 0; i < children.size(); i++){
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
		}

		//paralelizada 2
		void depth_first_searchP2(vector <int> &bestTour, vector <int> &tour, int totalCost){
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				#pragma omp for
				for(int i = 0; i < children.size(); i++){
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
		}
		//paralelizada 3
		void depth_first_searchP3(vector <int> &bestTour, vector <int> &tour, int totalCost){
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				#pragma omp for
				for(int i = 0; i < children.size(); i++){
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						#pragma omp critical
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
		}
		//paralelizada 4
		void depth_first_searchP4(vector <int> &bestTour, vector <int> &tour, int totalCost){
			#pragma omp single
			{
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				#pragma omp parallel for
				for(int i = 0; i < children.size(); i++){
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
			}
		}
		//paralelizada 5
		void depth_first_searchP5(vector <int> &bestTour, vector <int> &tour, int totalCost){
			if(tour.size() == n_cities){
				bestTour = tour;
			}else{
				#pragma omp for
				for(int i = 0; i < children.size(); i++){
					#pragma omp master
					if(feasible(tour, children[i]->data, totalCost)){
						tour.push_back(children[i]->data);
						totalCost += children[i]->cost;
						depth_first_search(bestTour, tour, totalCost);
						tour.pop_back();
					}
				}			
			
			}
		}
		//paralelizada 6
		void depth_first_search2P1(stack <int> S, vector <int> &tour){
			#pragma omp parallel for
			for(int i = n_cities-1; i >= 0; i--){
				S.push(i);
			}
			int totalCost = 0;
			#pragma omp parallel for
			while(!S.empty()){
				int city = S.top();
				if(city == NC){
					tour.pop_back();
				}
				else{
					tour.push_back(data);
					if(tour.size() == n_cities){
						if(feasible(tour, city, totalCost)){
							tour.push_back(data);
						}
						tour.pop_back();

					}else{
						S.push(NC);
						#pragma omp parallel for
						for(int i = n_cities-1; i >= 0; i++){
							if(feasible(tour, city, totalCost))
								S.push(city);
						}
					}
				}
			}
		}

		//paralelizada 7
		void depth_first_search2P2(stack <int> S, vector <int> &tour){
			#pragma omp for
			for(int i = n_cities-1; i >= 0; i--){
				S.push(i);
			}
			int totalCost = 0;
			#pragma omp for
			while(!S.empty()){
				int city = S.top();
				if(city == NC){
					tour.pop_back();
				}
				else{
					tour.push_back(data);
					if(tour.size() == n_cities){
						if(feasible(tour, city, totalCost)){
							tour.push_back(data);
						}
						tour.pop_back();

					}else{
						S.push(NC);
						#pragma omp for
						for(int i = n_cities-1; i >= 0; i++){
							if(feasible(tour, city, totalCost))
								S.push(city);
						}
					}
				}
			}
		}

		//paralelizada 8
		void depth_first_search2P3(stack <int> S, vector <int> &tour){
			#pragma omp for
			for(int i = n_cities-1; i >= 0; i--){
				S.push(i);
			}
			int totalCost = 0;
			#pragma omp for
			while(!S.empty()){
				int city = S.top();
				if(city == NC){
					tour.pop_back();
				}
				else{
					tour.push_back(data);
					if(tour.size() == n_cities){
						if(feasible(tour, city, totalCost)){
							tour.push_back(data);
						}
						tour.pop_back();

					}else{
						S.push(NC);
						#pragma omp critical
						for(int i = n_cities-1; i >= 0; i++){
							if(feasible(tour, city, totalCost))
								S.push(city);
						}
					}
				}
			}
		}
		//paralelizada 9
		void depth_first_search2P4(stack <int> S, vector <int> &tour){
			#pragma omp single
			for(int i = n_cities-1; i >= 0; i--){
				S.push(i);
			}
			int totalCost = 0;
			#pragma omp for
			while(!S.empty()){
				int city = S.top();
				if(city == NC){
					tour.pop_back();
				}
				else{
					tour.push_back(data);
					if(tour.size() == n_cities){
						if(feasible(tour, city, totalCost)){
							tour.push_back(data);
						}
						tour.pop_back();

					}else{
						S.push(NC);
						#pragma omp master
						for(int i = n_cities-1; i >= 0; i++){
							if(feasible(tour, city, totalCost))
								S.push(city);
						}
						#pragma omp barrier
					}
				}
			}
		}
};	

int _tmain(int argc, _TCHAR* argv[])
{
	srand((unsigned)time(NULL));
	Matrix *digraph = new Matrix(n_cities, n_cities);
	digraph->print();
	Node *tree = new Node();
	tree->buildTree(digraph, hometown, 0);
	tree->print();
	vector <int> tour; //almacena el tour en el depth search
	vector <int> bestTour;
	int totalCost = 0;
	stack <int> S;
	tree->depth_first_searchP4(bestTour, tour, totalCost);
	tree->depth_first_search2P3(S, tour);
	return 0;
}

