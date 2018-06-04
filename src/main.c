  #include <stdio.h>
  #include <stdlib.h>
  #include <mpi.h>
  #include <math.h>

  #define SOFTENING 1e-9f


  typedef struct{
      float x;
      float y;
      float z;
      float vx;
      float vy;
      float vz;
  }Body; 

  MPI_Datatype BodyMPI;

  Body *bodies;

  //MPI ENV
  int  my_rank; /* rank of process */
  int  p;       /* number of processes */
  MPI_Status status;   /* return status for receive */


  //prototype
  void checkArgs(int argc, char *argv[]);
  void initBodies();
  void printBodies(Body *body,int lenght);
  void bodyForce(Body *bodyPart,int lenght,int start);
  void updatePositions(Body *bodyPart, int lenght,int start);
  //args
  int particelle; //numero di particelle da computare, inserito dall'utente
  int numIter; //iterazioni da fare, inserito dall'utente

  const double dt=0.1f; //tempo dell'iterazione

  int main(int argc, char* argv[]){


    /* start up MPI */
    MPI_Init(&argc, &argv);
    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    checkArgs(argc,argv);
    
    bodies = malloc(particelle * sizeof(Body));//alloca l'insieme delle particelle
    MPI_Type_contiguous(6, MPI_FLOAT, &BodyMPI);//creiamo un tipo contiguoi in modo da porterlo utilizzare nelle comunicazioni MPI
    MPI_Type_commit(&BodyMPI);//rende effettiva la allocazione di memoria
    
    
    initBodies();//tutti i processori inizializzano il proprio bodies
    
    double tstart = MPI_Wtime();//tempo iniziale
    int reminder = particelle % p;
    int partition = particelle / p;
    int p1= partition + 1;
    
    if(reminder==0){
      //se il reminder è pari a 0 inviamo a tutti la stessa porzione di array sulla quale bisogna chiamare bodyForce
      MPI_Scatter(bodies,partition,BodyMPI,bodies,partition,BodyMPI,0,MPI_COMM_WORLD);
      //la funzione bodyForce calcola lo spostamento delle particelle su un body grande partition partendo da partition * myrank
      bodyForce(bodies,partition,partition*my_rank);
      
    }else{
      if(my_rank==0){//MASTER
        int k=1;
        int start=p1; //da quale posizione deve compinciare la computazione di bodyForce
        //invia a tutti i processori, tranne al processore 0, la parte di bodies da computare
        for(int i=1;i<p;i++){
          if(k<reminder){
            //se il resto è diverso da 0, allora fino a quando k sarà minore di reminder a ogni processore viene inviato un bodies grande p1
            MPI_Send(&bodies[start],p1,BodyMPI,i,0,MPI_COMM_WORLD);
            k++;
            start+=p1;//aggiorno il punto di partenza
          }else{
         //una volta che k è diventato uguale a reminder vuol dire che non abbiamo più resto quindi possiamo inviare partition
          MPI_Send(&bodies[start],partition,BodyMPI,i,0,MPI_COMM_WORLD);
         
            start+=partition;
          }
        }
     //il MASTER calcola gli spostamenti delle particelle, partendo da 0
     bodyForce(bodies,p1,0);

      }else{//SLAVE
        //se my_rank è minore di reminder vuol dire che devo ricevere un bodies grande p1
        if(my_rank<reminder){
         MPI_Recv(bodies,p1,BodyMPI,0,0,MPI_COMM_WORLD,&status);
         bodyForce(bodies,p1,p1*my_rank);
          
        }
    else{
      //altrimenti vuol dire che devo ricevere un bodies grande partition
      MPI_Recv(bodies,partition,BodyMPI,0,0,MPI_COMM_WORLD,&status);
      bodyForce(bodies,partition,partition*my_rank);
      }//fine SLAVE
    }//else reminder 0
  }

    MPI_Barrier(MPI_COMM_WORLD);//aspetto che tutti i processori abbiano terminato la computazione
    //il master stampa i tempi finali e bodies
    if(my_rank==0){
     double tend = MPI_Wtime();// tempo finale della computazione della nbody simulation
     printBodies(bodies,particelle);//stampa bodies
     double totTime = tend-tstart;// tempo totale
     printf("totTime: %f\n",totTime);
    }
    
    MPI_Type_free(&BodyMPI);// libera la memoria allocata con MPI_Type_contiguous
    free(bodies);//dealloca bodies
    MPI_Finalize();
    return 0;
  }

  /**
    *cotrolla gli argomenti passati da riga di comando
    */
  void checkArgs(int argc, char *argv[]){
    if(argc != 3){
      printf("Inserisci parametri\n");
      printf("Numero di particelle da genereare\n");
      printf("Numero di iterazioni da fare\n");
    }else{
          particelle = atoi(argv[1]);
          numIter = atoi(argv[2]);
          if(particelle<=0 || numIter<=0){
            MPI_Finalize();
            exit(0);
          }
    }
  }

  /**
    *inizializza bodies
    */
  void initBodies(){
    for(int i=0; i< particelle; i++){
      bodies[i].x=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
      bodies[i].y=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
      bodies[i].z=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
      bodies[i].vx=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
      bodies[i].vy=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
      bodies[i].vz=2.0f * (rand() / (float) RAND_MAX)  - 1.0f;
    }

  }
  /**
    *calcola gli spostamenti delle particelle e aggiorna le posizioni
    */
  void bodyForce(Body *bodyPart, int lenght, int start) {
    printf("Sono %d e ricevo %d particella e parto da %d \n",my_rank,lenght,start );
    //calcola gli spostamenti sulla porzione di bodies inviata da un processore
    for (int i = start; i < start + lenght; i++) { 
      float Fx = 0.0f;
      float Fy = 0.0f;
      float Fz = 0.0f;
      for(int it= 0; it < numIter;it++){
      for (int j = 0; j < particelle; j++) {
        float dx = bodyPart[j].x - bodyPart[i].x;
        float dy = bodyPart[j].y - bodyPart[i].y;
        float dz = bodyPart[j].z - bodyPart[i].z;
        float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * invDist3;
        Fy += dy * invDist3;
        Fz += dz * invDist3;
      }

      bodyPart[i].vx += dt*Fx;
      bodyPart[i].vy += dt*Fy;
      bodyPart[i].vz += dt*Fz;
      
    } 
    for (int i = start ; i < start + lenght; i++) { 
        bodyPart[i].x += bodyPart[i].vx*dt;
        bodyPart[i].y += bodyPart[i].vy*dt;
        bodyPart[i].z += bodyPart[i].vz*dt;
      }
    }
    
  }

  /**
    * stampa bodies
    */
  void printBodies(Body *body, int lenght){
        for(int i=0;i<lenght;i++){
        printf("--------------------------------------------------%d--------------------------------------------------\n",i);
        printf("x= %f\ty= %f\tz= %f\tvx= %f\tvy= %f\tvz= %f\t\n",body[i].x,body[i].y,body[i].z,body[i].vx,body[i].vy,body[i].vz);
        printf("--------------------------------------------------%d--------------------------------------------------\n\n",i);
      }
  }