#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <time.h>
#include <vector>
#include <math.h>

using namespace std;

#define RANDOM_LIMIT 9
#define n_cities 5
#define hometown 0 

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
};

int _tmain(int argc, _TCHAR* argv[])
{
	srand((unsigned)time(NULL));
	Matrix *digraph = new Matrix(n_cities, n_cities);
	digraph->print();
	Node *tree = new Node();
	tree->buildTree(digraph, hometown, 0);
	tree->print();
	
	
	return 0;
}

