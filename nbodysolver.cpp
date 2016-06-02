#include <iostream>
#include <omp.h>
#include <time.h>
#include <vector>
#include <math.h>

using namespace std;
#define RANDOM_LIMIT 100000
#define G 6.673e-11 //6.673 * pow(10.0, -11.0)

class Position{
	double x;
	double y;
	public:
		Position(){
			x = (double)(rand()%RANDOM_LIMIT);
			y = (double)(rand()%RANDOM_LIMIT);
		}
		void setX(int x1){
			x = x1;
		}
		void setY(int y1){
			y = y1;
		}
		double getX(){
			return x;
		}
		double getY(){
			return y;
		}
		void addToX(double p_x){
			x += p_x;
		}
		void addToY(double p_y){
			y += p_y;
		}
};

class Particle{
	public:
		Position *position;
		double force[2];
		double mass;
		double velocity[2];
		double acceleration;
	
		Particle(){
			position = new Position();
			mass = 10.0; //mass = (double)(rand()%10);
		}
		void setForce(double f_x, double f_y){
			force[0] = f_x;
			force[1] = f_y;
		}
		void setMass(double m){
			mass = m;
		}
		Position* getPosition(){
			return position;
		}
		double* getForce(){
			return force;
		}
		double getMass(){
			return mass;
		}
		void addToForce(double f_x, double f_y){
			force[0] += f_x;
			force[1] += f_y;
		}
		void minusToForce(double f_x, double f_y){
			force[0] -= f_x;
			force[1] -= f_y;
		}
};

//como se especifica el usuario ingresa los siguientes valores
static long num_particles = 1000;
double delta_t = 10.0;
int timestep = 10;

//version secuencial 1
void nbodySolverSeq(vector <Particle*> &particles){ 
	int i, j;
	double x_diff, y_diff;
	double dist, dist_cubed;
	for(i = 0; i < num_particles; i++){
		for(j = 0; j < num_particles; j++){
			if(i != j){ 
				x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
				y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
				dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
				dist_cubed = dist*dist*dist;
				double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
				double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

				particles[i]->minusToForce(force_ij_X, force_ij_Y);
				//particles[j]->minusToForce(force_ij_X, force_ij_Y);
			}
		}
		
	}
	for(i = 0; i < num_particles; i++){
		particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
		particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
		particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
		particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
	}
}
//version secuencial 2 reducida
void nbodySolverSeq2(vector <Particle*> &particles){ 
	int i, j;
	double x_diff, y_diff;
	double dist, dist_cubed;
	for(i = 0; i < num_particles; i++){
		for(j = 0; j < i; j++){
				x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
				y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
				dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
				dist_cubed = dist*dist*dist;
				double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
				double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

				particles[i]->addToForce(force_ij_X, force_ij_Y);
				particles[j]->minusToForce(force_ij_X, force_ij_Y);
				//particles[j]->minusToForce(force_ij_X, force_ij_Y);
		}		
	}
	for(i = 0; i < num_particles; i++){
		particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
		particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
		particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
		particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
	}
	/*cout<<"particle "<<j<<" has velocity = ("<<particles[0]->velocity[0]<<", "<<particles[0]->velocity[1]<<
					") and position ("<<particles[0]->position->getX()<<", "<<particles[0]->position->getY()<<")"<<endl;*/
}

//version paralela 1
void nbodySolverPar1(vector <Particle*> &particles){ 
    #pragma omp parallel for
	for(int i = 0; i < num_particles; i++){
		for(int j = 0; j < i; j++){
				double x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
				double y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
				double dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
				double dist_cubed = dist*dist*dist;
				double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
				double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

				particles[i]->addToForce(force_ij_X, force_ij_Y);
				particles[j]->minusToForce(force_ij_X, force_ij_Y);
				//particles[j]->minusToForce(force_ij_X, force_ij_Y);
		}		
	}
	#pragma omp parallel for
	for(int i = 0; i < num_particles; i++){
		particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
		particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
		particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
		particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
	}
}

//version paralela 1.1
void nbodySolverPar1_1(vector <Particle*> &particles){ 
    #pragma omp for
	for(int i = 0; i < num_particles; i++){
		for(int j = 0; j < i; j++){
				double x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
				double y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
				double dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
				double dist_cubed = dist*dist*dist;
				double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
				double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

				particles[i]->addToForce(force_ij_X, force_ij_Y);
				particles[j]->minusToForce(force_ij_X, force_ij_Y);
				//particles[j]->minusToForce(force_ij_X, force_ij_Y);
		}		
	}
	#pragma omp for
	for(int i = 0; i < num_particles; i++){
		particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
		particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
		particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
		particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
	}
}

//version paralela 2 evitando varias llamadas a particles
void nbodySolverPar2(vector <Particle*> &particles){ 
    #pragma omp parallel for
	for(int i = 0; i < num_particles; i++){
		for(int j = 0; j < i; j++){
				Particle *p_i = particles[i];
				Particle *p_j = particles[j];
				Position *pos_i = p_i->getPosition();
				Position *pos_j = p_j->getPosition();
				double m_i = p_i->getMass();
				double m_j = p_j->getMass();
				double x_diff = pos_i->getX() - pos_j->getX();
				double y_diff = pos_i->getY() - pos_j->getY();
				double dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
				double dist_cubed = dist*dist*dist;
				double force_ij_X = (G*m_i*m_j)/(dist_cubed*x_diff);
				double force_ij_Y = (G*m_i*m_j)/(dist_cubed*y_diff);

				particles[i]->addToForce(force_ij_X, force_ij_Y);
				particles[j]->minusToForce(force_ij_X, force_ij_Y);
				//particles[j]->minusToForce(force_ij_X, force_ij_Y);
		}		
	}
	#pragma omp parallel for
	for(int i = 0; i < num_particles; i++){ 
		particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
		particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
		particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
		particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
	}
}
//versión paralela 3 paralelizando sólo el primer loop de timesteps y usando la primera versión secuencial
void nBodySolverPar3(vector <Particle*> &particles){
	#pragma omp parallel for
	for(int i = 0; i < timestep; i++){
		nbodySolverSeq(particles);
	}
}

//versión paralela 4 paralelizando todos los loops
void nBodySolverPar4(vector <Particle*> &particles){
	#pragma omp parallel for
	for(int i = 0; i < timestep; i++){
		nbodySolverPar1(particles);
	}
}

//version paralela 5 usando #pragma omp parallel en el 1er loop y  #pragma omp for en los demás 
void nBodySolverPar5(vector <Particle*> &particles){
	#pragma omp parallel
	for(int i = 0; i < timestep; i++){
		nbodySolverPar1_1(particles);
	}
}

//version paralela 6 usando versión reducida con un loop adicional para inicializar las fuerzas
void nBodySolverPar6(vector <Particle*> &particles){
	#pragma omp parallel
	for(int i = 0; i < timestep; i++){
		#pragma omp for
		for(int i = 0; i < num_particles; i++){
			particles[i]->setForce(0.0,0.0);
		}
		#pragma omp for
		for(int i = 0; i < num_particles; i++){
			for(int j = 0; j < i; j++){
					double x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
					double y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
					double dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
					double dist_cubed = dist*dist*dist;
					double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
					double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

					particles[i]->addToForce(force_ij_X, force_ij_Y);
					particles[j]->minusToForce(force_ij_X, force_ij_Y);
					//particles[j]->minusToForce(force_ij_X, force_ij_Y);
			}		
		}
		#pragma omp for
		for(int i = 0; i < num_particles; i++){
			particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
			particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
			particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
			particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
		}
	}
}

//version paralela 7 usando #pragma omp critical
void nBodySolverPar7(vector <Particle*> &particles){
	#pragma omp parallel
	for(int i = 0; i < timestep; i++){
		#pragma omp for
		for(int i = 0; i < num_particles; i++){
			for(int j = 0; j < i; j++){
					double x_diff = particles[i]->getPosition()->getX() - particles[j]->getPosition()->getX();
					double y_diff = particles[i]->getPosition()->getY() - particles[j]->getPosition()->getY();
					double dist = sqrt((double)(x_diff*x_diff + y_diff*y_diff));
					double dist_cubed = dist*dist*dist;
					double force_ij_X = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*x_diff);
					double force_ij_Y = (G*particles[i]->getMass()*particles[j]->getMass())/(dist_cubed*y_diff);

					particles[i]->addToForce(force_ij_X, force_ij_Y);
					particles[j]->minusToForce(force_ij_X, force_ij_Y);
					//particles[j]->minusToForce(force_ij_X, force_ij_Y);
			}		
		}
		#pragma omp critical
		for(int i = 0; i < num_particles; i++){
			particles[i]->position->addToX(delta_t*particles[i]->velocity[0]);
			particles[i]->position->addToY(delta_t*particles[i]->velocity[1]);
			particles[i]->velocity[0] += delta_t/(particles[i]->mass*particles[i]->force[0]);
			particles[i]->velocity[1] += delta_t/(particles[i]->mass*particles[i]->force[1]);
		}
	}
}

int main(){
	srand((unsigned)time(NULL));
	vector <Particle*> particles(num_particles);
	for(int i = 0; i < num_particles; i++){ //inicializacion de las fuerzas de cada particula con valor 0.0
		particles[i] = new Particle(); //al crear la particula inmediatamente el constructor setea valores aleatorios en las pocisiones (x,y) del plano
		particles[i]->setForce(0.0, 0.0); //fuerza en x, y
		particles[i]->velocity[0] = 0.1; //velocidad inicial en x
		particles[i]->velocity[1] = 0.1; //velocidad inicial en y
	}
	clock_t start = clock();
	/*
	for(int i = 0; i < timestep; i++){
		//for(int j = 0; j < 1; j++){
		//	cout<<"particle "<<j<<" has velocity = ("<<particles[j]->velocity[0]<<", "<<particles[j]->velocity[1]<<
		//		  ") and position ("<<particles[j]->position->getX()<<", "<<particles[j]->position->getY()<<")"<<endl;
		//}
		nbodySolverPar1(particles);
	}*/
	nBodySolverPar7(particles);

	printf("Tiempo transcurrido: %f", ((double)clock() - start) / CLOCKS_PER_SEC);

	return 0;
}

//double omp_get_wtime
/*
static long num_steps = 1000000000;
int main(){ 		
	clock_t start = clock();
	int i;
	double x, pi, sum = 0.0;
	double step = 1.0/(double) num_steps;
  
	for (i=0; i < num_steps; ++i) {
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}  
	pi = step * sum;
	cout<<pi<<endl;
	printf("Tiempo transcurrido: %f", ((double)clock() - start) / CLOCKS_PER_SEC);
}
*/

//version paralela
/*
static long num_steps = 1000000000;
#define NUM_THREADS 4

long v1[] = {0, 250000000, 500000000, 750000000};
long v2[] = {250000000, 500000000, 750000000, 1000000000};

int main(){
	clock_t start = clock();
	double step = 1.0/(double)num_steps;
	double sum[NUM_THREADS];
	double pi = 0.0;

	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int id = omp_get_thread_num();
		int i;
		double x;
		int nts = omp_get_num_threads();
		//cout<<"ss "<<nts<<endl;
		double sum_tmp = 0.0;
		for(i = id; i < num_steps; i+=NUM_THREADS){
			x = (i+0.5)*step;
			sum_tmp += 4.0/(1.0+x*x);
		}
		sum[id] = sum_tmp;
		//cout<<"acabo el for con valor "<<i<<" del id "<<id<<" r = "<<sum[id]<<endl;
	}
	double s = 0.0;
	for(int i = 0; i < NUM_THREADS; i++)
		pi += step*sum[i];
	//pi = step * s;
	cout<<pi<<endl;
	printf("Tiempo transcurrido: %f", ((double)clock() - start) / CLOCKS_PER_SEC);
}
*/
/*
//version paralela usando parallel for
static long num_steps = 1000000000;
#define NUM_THREADS 4
int main(){ 		
	clock_t start = clock();
	int i;
	double x, pi, sum = 0.0;
	double step = 1.0/(double) num_steps;
	#pragma omp parallel for
	for (i=0; i < num_steps; ++i) {
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}  
	
	pi = step * sum;
	cout<<pi<<endl;
	printf("Tiempo transcurrido: %f", ((double)clock() - start) / CLOCKS_PER_SEC);
}
*/
/*
int main(){
	double A[1000];
	//omp_set_num_threads(4);
	#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		cout<<"mi id es "<<id<<endl;
	}

	return 0;
}
*/

/*
int main(){
	#pragma omp parallel
	{
	int id = omp_get_thread_num();
	cout<<"Hello "<<id<<endl;
	cout<<"World "<<id<<endl;
	}
	return 0;
}*/