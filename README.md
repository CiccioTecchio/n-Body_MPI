# N-body Simulation

***
## Programmazione Concorrente, Parallela e su Cloud

### Università  degli Studi di Salerno

#### *Anno Accademico 2017/2018*

**Professore:** *Vittorio Scarano,*
 **Dottor:** *Carmine Spagnuolo,*
 **Studente:** *Francesco Vicidomini*


---

## Problem Statement  

 Nel n-body problem, abbiamo bisogno di trovare le posizioni e le velocità di una collezione di particelle che interagistrono fra loro per un determinato periodo di tempo.
 Per esempio, un astrofisico è interessato a conoscere la posizione e alla velocità di una collezione di stelle. Un n-body solver è un programma che cerca la soluzione a un n-body problem simulando il movimento delle particelle. All'n-body solver vengono dati in input il numero di particelle, in maniera casuale verranno assegnate le posizioni nelle spazio e le relative velocità e il numero di iteazioni che deve simulare. L'output sarà la posizione e la velocità di ogni particella alla fine di un determinato numero di iterazioni specificato dall'utente.

## Soluzione proposta

La soluzione proposta considera solo l'approccio n^2 rispetto al numero di particelle e al numero di iterazioni dato dall'utente.
La comunicazione è avvenuta usando sia la comunicazione collettiva con la funzione **MPI_Scatter** che la comunicazione Point to Point con le funzioni **send** e **recv**
I test sono stati effettuati sulle istanze di AWS **m4.large**. 

### Implementazione

L'obiettivo del lavoro svolto è stato quello di parallelizzare l'algoritmo dell' n-body simulationio, partizionando equamente il lavoro tra i processi coinvolti.
L'approccio utilizzato al fine di ottenere una distribuzione quanto più equa possibile è quello riportato di seguito.
#### Descrizione variabili
**Body:** struttura utilizzata per rappresentare una singola particella
```c
  typedef struct{
      float x;
      float y;
      float z;
      float vx;
      float vy;
      float vz;
  }Body; 
```
**Variabili MPI:**
```c
int my_rank; //rank del processore
int p; //numero totale di processori
MPI_Datatype BodyMPI; //dichiariamo il tipo di dato BodyMPI 
```
**Variabili passate da riga di comando**
```c
int particelle; //numero di particelle da creare
int numIter //numero di iterazioni da eseguire
```
```c
const double dt=0.1f; //tempo dell'iterazione
```
#### checkArgs
Controlla se gli argomenti passati da riga di comando sono corretti.
```c
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
```
#### initBodies
Tutti i processori inizializzano il loro *bodies*
```c
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
  ```

#### bodyForce: ####
calcola gli spostamenti sulla porzione di bodies inviata da un processore.
 - bodyPart: rappresenta la porzione di bodies che un singolo processore deve computare;
 - lenght: quanto è grande bodyPart;
 - start: da quale indice di bodies deve iniziare a computare.

```c
  void bodyForce(Body *bodyPart, int lenght, int start) {
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
```
#### printBodies: ####
Stampa tutte le componenti di un Body grande lenght
```c
  void printBodies(Body *body, int lenght){
        for(int i=0;i<lenght;i++){
        printf("--------------------------------------------------%d--------------------------------------------------\n",i);
        printf("x= %f\ty= %f\tz= %f\tvx= %f\tvy= %f\tvz= %f\t\n",body[i].x,body[i].y,body[i].z,body[i].vx,body[i].vy,body[i].vz);
        printf("--------------------------------------------------%d--------------------------------------------------\n\n",i);
      }
  }
```

### Il main:

La prima operazione effettuata nel main è quella di inizializzare l'ambiente MPI e di controllare se i parametri passati dell'utente sono corretti.
```c
 /* start up MPI */
    MPI_Init(&argc, &argv);
    /* find out process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    checkArgs(argc,argv);
```

Fatto ciò viene allocata la memoria necessaria per contenere un insieme di particelle.
```c
bodies = malloc(particelle * sizeof(Body));
```
e con la MPI_Type_contiguous rendiamo le variabili di tipo BodyMPI accessibili in maniera contigua
```c
MPI_Type_contiguous(6, MPI_FLOAT, &BodyMPI);
MPI_Type_commit(&BodyMPI);
```

Inizializziamo bodies
```c
initBodies()
```
Inizializziamo il tempo di esecuzione, il reminder, la partition e p1.
```c
    double tstart = MPI_Wtime();
    int reminder = particelle % p;
    int partition = particelle / p;
    int p1= partition + 1;

```
Se il reminder è pari a 0 allora possiamo effettuare una chiamata collettiva perchè il master e gli slave devono ricevere tutti la stessa parzione di bodies che è partition.
```c
MPI_Scatter(bodies,partition,BodyMPI,bodies,partition,BodyMPI,0,MPI_COMM_WORLD);
bodyForce(bodies,partition,partition*my_rank);
```
#### L'invio
Se il reminder è diverso da 0 dobbiamo utilizzare le chiamate Point 2 Point perchè la grandezza del bodies non è sempre fissa dato che dobbiamo inviare delle parti grandi p1.
Il master inizializza k e start, fatto ciò controlla se k è minore di reminder inviara una porzione di bodies grande p1 al processore i, incrementa k e aggiorna start; altrimenti se k è maggiore o uguale a reminder vuol dire che non dobbiamo più inviare parti grandi p1 perchè il resto è esaurito quindi inviamo una parte di bodies grande partition e aggiorniamo start. In fine il master chiama il suo bodyForce partendo da 0

```c
if(my_rank==0){//MASTER
  int k=1;
  int start=p1;
  for(int i=1;i<p;i++){
 if(k<reminder){
         MPI_Send(&bodies[start],p1,BodyMPI,i,0,MPI_COMM_WORLD);
         k++;
         start+=p1;
 }else{
MPI_Send(&bodies[start],partition,BodyMPI,i,0,MPI_COMM_WORLD);
start+=partition;
          }
        }
     bodyForce(bodies,p1,0);
```
#### La ricezione
Alcuni Slave dovranno ricevere un bodies grande p1 altri dovranno ricevere un bodies grande partition.
Se my_rank è minore di reminder vuol dire che devo ricevere una porzione grande p1 e sulla porzione ricevuta chiamerò bodyForce. Altrimenti devo ricevere un bodies grande partition e chiamare bodyForce
```c
else{//SLAVE
        if(my_rank<reminder){
         MPI_Recv(bodies,p1,BodyMPI,0,0,MPI_COMM_WORLD,&status);
         bodyForce(bodies,p1,p1*my_rank);   
        }
    else{
MPI_Recv(bodies,partition,BodyMPI,0,0,MPI_COMM_WORLD,&status);
bodyForce(bodies,partition,partition*my_rank);
      }//fine SLAVE
```
#### Conclusione
Usiamo una MPI_Barrier per assicurarci che tutti i processori abbiano completato la computazione, il MASTER stampa bodies eil tempo totale della computazione.
```c
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank==0){
     double tend = MPI_Wtime();
     printBodies(bodies,particelle);
     double totTime = tend-tstart;
     printf("totTime: %f\n",totTime);
    }
    MPI_Type_free(&BodyMPI);
    free(bodies);
    MPI_Finalize();
    return 0;
```
### Testing

I test sono stati effettuati sulle istanze __m4.large__ (2 core) di Amazon Web Services.  
Durante la fase di testing si è tenuto conto sia di strong scaling che di weak scaling

Risorse massime utilizzate:

* 8 Istanze EC2 m2.xlarge **StarCluster-Ubuntu-12.04-x86_64-hvm** - _ami-52a0c53b_
* 16 processori (2 core per Istanza).
## Strong Scaling

Nella fase di testing che ha tenuto in considerazione lo strong scaling sono state utilizzate 50.000 particelle e 20 iterazioni.
Nello strong scaling infatti il numero di particelle resta invariato, quello che cambia è il numero di processori.
Nella figura in basso è possibile osservare i risultati di questa fase di testing.

![enter image description here](https://lh3.googleusercontent.com/WseuRAGXmcsSMqR1uOoBBgafdf054H4gFtPFai010UmNO7RQlHMR_JrhdY2aeBzAFaWHWFC7mAs "Strong test")


## Weak Scaling

La fase di testing che ha tenuto in considerazione il weak scaling è stata svolta in due parti.
Inizialmente sono state utilizzate 3000 particelle e 20 iterazione per processo.
In seguito 10000 particelle e 20 iterazioni per processo.
Nel weak scaling infatti il numero di particelle cresce in maniera proporzionale al numero di processori.
Nella figura in basso è possibile osservare i risultati di questa fase di testing.

![](https://lh3.googleusercontent.com/ypY7MTn4XMy9GyO5pMpu23WS_S89A6PiDUu27z-nRKp6kYZbLkC7aDTGmraprkRlmzD35hhNvQY "Weak")
