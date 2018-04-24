#include "define.h"

#define MASTER 0
#define DIM 2
#define X 0
#define Y 1
#define FILENAME "input_file_8.txt"
#define DEBUG_GET_ARGUMENT 1
//#define DEBUG_OUTPUT_STATE 1
//#define DEBUG_FORCES_BEFORE 0
//#define DEBUG_FORCES_AFTER 0
//#define DEBUG_READ_FILE 0
//#define DEBUG_UPDATE_BEFORE 0
//#define DEBUG_UPDATE_AFTER 0
#define GENERATE_INPUT_FILE 1
#define GENERATE_OUTPUT_FILE 1

typedef double vector[DIM];
const double G = 6.673e-11;

int my_rank;
int num_processes; // number of processes

MPI_Datatype vectorMPI;

vector *velocities = NULL;

int img_width = 2550;
int img_height = 1440;

void Get_Input_Arguments(int argc, char* argv[], int* numParticlesLight, int* numParticleMedium, int* numParticleHeavy, int* num_particles, int* num_steps, double* delta_t, int* output_freq);
void Read_File_Init_Conditions(char* filaName, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void Generate_Init_Conditions(int types[], double masses[], vector positions[], vector my_velocities[], int numParticlesLight, int numParticleMedium, int numParticleHeavy, int num_particles, int chunk);
void Output_State(double time, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void update_state(unsigned char* image, int types[], double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);
void Compute_Force(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk);
void Update_Particles(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double delta_t);
void Generate_Output_File(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk);

//main program
int main(int argc, char* argv[]){
    srand (1234);

    int numParticlesLight, numParticleMedium, numParticleHeavy, num_particles;

    int chunk, num_steps, steps, subSteps, my_particles, output_freq, *types;
    double delta_t, time, *masses;

    vector* my_positions, * positions, * my_velocities, * my_forces;

    double start_time, end_time;

    // MPI bookkeeping
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    // Get input
    //if (argc != 10) {
    //    printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
    //}
    if(my_rank == MASTER){
        numParticlesLight = atoi(argv[1]);
        numParticleMedium = atoi(argv[2]);
        numParticleHeavy = atoi(argv[3]);
        num_steps = atoi(argv[4]);
        subSteps = atoi(argv[5]); // Doesn't use right now
        delta_t = strtod(argv[6], NULL);
        //img_width = atoi(argv[7]);
        //img_height = atoi(argv[8]);

        num_particles = numParticlesLight + numParticleMedium + numParticleHeavy;
        output_freq = 1; // default
    }

    MPI_Bcast(&num_particles, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&num_steps, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&delta_t, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&output_freq, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    #  ifdef DEBUG_GET_ARGUMENT
    if (my_rank == 0) {
        printf("num_particles = %d\n", num_particles);
        printf("num_steps = %d\n", num_steps);
        printf("delta_t = %e\n", delta_t);
        printf("output_freq = %d\n", output_freq);
    }
    #  endif

    // Distribute
    int have_extra = num_particles % num_processes;

    int padding_num = (have_extra == 0) ? 0 : num_processes - have_extra;
    num_particles += padding_num;
    chunk = num_particles / num_processes;


    types = (int *)malloc(num_particles*sizeof(int));
    masses = (double *)malloc(num_particles*sizeof(double));
    positions = (vector *)malloc(num_particles*sizeof(vector));
    my_forces = (vector *)malloc(chunk*sizeof(vector));
    my_positions = positions + my_rank*chunk;
    my_velocities = (vector *)malloc(chunk*sizeof(vector));

    if (my_rank == MASTER) {
        velocities = (vector *)malloc(num_particles*sizeof(vector));
    }

    MPI_Type_contiguous(DIM, MPI_DOUBLE, &vectorMPI);
    MPI_Type_commit(&vectorMPI);

    // Generate particles
    Generate_Init_Conditions(types, masses, positions, my_velocities, numParticlesLight, numParticleMedium, numParticleHeavy, num_particles, chunk);

    start_time = MPI_Wtime();

    #ifdef DEBUG_OUTPUT_STATE
    Output_State(0.0, masses, positions, my_velocities, num_particles, chunk);
    #endif

    unsigned char* image = (unsigned char*) malloc(img_width * img_height * 3 * sizeof(unsigned char));
    for (steps = 1; steps <= num_steps; steps++) {
        time = steps*delta_t;
        for (my_particles = 0; my_particles < chunk; my_particles++) {
            if (types[my_particles] == DUMMY) {
                continue;
            }
            Compute_Force(my_particles, masses, my_forces, positions, num_particles, chunk);
        }

        for (my_particles = 0; my_particles < chunk; my_particles++) {
            if (types[my_particles] == DUMMY) {
                continue;
            }
            Update_Particles(my_particles, masses, my_forces, my_positions, my_velocities, num_particles, chunk, delta_t);
        }

        MPI_Allgather(MPI_IN_PLACE, chunk, vectorMPI, positions, chunk, vectorMPI, MPI_COMM_WORLD);

        // Save BMP
        initilize_img(image, img_width, img_height);
        update_state(image, types, masses, positions, my_velocities, num_particles, chunk);

        // Save the image
        char img_name[100];
        char img_full_path[100];
        snprintf(img_name, sizeof(char) * 32, "_%05i.bmp", steps);

        strcpy(img_full_path, "../images/out");
        strcat(img_full_path, img_name);
        saveBMP(img_full_path, image, img_width, img_height);

        #ifdef DEBUG_OUTPUT_STATE
        if (steps % output_freq == 0)
            Output_State(time, masses, positions, my_velocities, num_particles, chunk);
        #endif

    }

    #ifdef GENERATE_OUTPUT_FILE
    Generate_Output_File(masses, positions,  my_velocities,  num_particles, chunk);
    #endif

    end_time=MPI_Wtime();
    if(my_rank==MASTER){
        printf("Time passed = %.5f seconds\n", end_time-start_time);
        // Output to csv. num_particles, num_steps, num_processors, time taken (s)
        fprintf(stderr, "%d,%d,%d,%.5f\n",num_particles,num_steps,num_processes,end_time-start_time);
    }

    MPI_Type_free(&vectorMPI);
    free(masses);
    free(positions);
    free(my_forces);
    free(my_velocities);
    if (my_rank == MASTER)
        free(velocities);

    MPI_Finalize();
    return 0;
}


void Generate_Init_Conditions(int types[], double masses[], vector positions[], vector my_velocities[], int numParticlesLight, int numParticleMedium, int numParticleHeavy, int num_particles, int chunk){
    int part;
    int line;
    //double mass = 5.0e24;
    double mass = 5;
    //double gap = 1.0e5;
    double gap = 10;

    if (my_rank == MASTER) {
        // light particles
        //for (part = 0; part < num_particles; part++) { // BUG ?
        for (part = 0; part < numParticlesLight; ++part) {
            types[part] = LIGHT;
            masses[part] = massLightMin + (massLightMax - massLightMin) * drand48();
            positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
            positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

            // generate velocities
            double max_v, min_v, rand_v, rand_deg;
            max_v = velocityLightMin;
            min_v = velocityLightMax;
            rand_v = min_v + (max_v - min_v) * drand48();
            rand_deg = TWO_PI * drand48();
            velocities[part][X] = rand_v * cos(rand_deg);
            velocities[part][Y] = rand_v * sin(rand_deg);
        }

        // medium particles
        for (part = numParticlesLight; part < numParticlesLight+numParticleMedium; ++part) {
            types[part] = MEDIUM;
            masses[part] = massMediumMin + (massMediumMax - massMediumMin) * drand48();
            positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
            positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

            // generate velocities
            double max_v, min_v, rand_v, rand_deg;
            max_v = velocityMediumMin;
            min_v = velocityMediumMax;
            rand_v = min_v + (max_v - min_v) * drand48();
            rand_deg = TWO_PI * drand48();
            velocities[part][X] = rand_v * cos(rand_deg);
            velocities[part][Y] = rand_v * sin(rand_deg);
        }

        // heavy particles
        for (part = numParticlesLight+numParticleMedium; part < numParticlesLight+numParticleMedium+numParticleHeavy; ++part) {
            types[part] = HEAVY;
            masses[part] = massHeavyMin + (massLightMax - massHeavyMin) * drand48();
            positions[part][X] = POS_MIN_X + (POS_MAX_X - POS_MIN_X) * drand48();
            positions[part][Y] = POS_MIN_Y + (POS_MAX_Y - POS_MIN_Y) * drand48();

            // generate velocities
            double max_v, min_v, rand_v, rand_deg;
            max_v = velocityHeavyMin;
            min_v = velocityHeavyMax;
            rand_v = min_v + (max_v - min_v) * drand48();
            rand_deg = TWO_PI * drand48();
            velocities[part][X] = rand_v * cos(rand_deg);
            velocities[part][Y] = rand_v * sin(rand_deg);
        }

        // padding particles
        for (part = numParticlesLight+numParticleMedium+numParticleHeavy; part < num_particles; ++part) {
            types[part] = DUMMY;
            masses[part] = 0;
            positions[part][X] = 0;
            positions[part][Y] = 0;

            // generate velocities
            velocities[part][X] = 0;
            velocities[part][Y] = 0;
        }
        
        #ifdef GENERATE_INPUT_FILE
        FILE *fp_generate = fopen("generate_input.txt","w+");
        if(!fp_generate){
            printf("Cannot create input file <%s>\n", "generate_input");
            exit(1);
        }

        for(line=0; line<num_particles; line++){
            fprintf(fp_generate, "%lf %lf %lf %lf %lf\n", masses[line], positions[line][X], positions[line][Y], velocities[line][X], velocities[line][Y]);
        }

        fclose(fp_generate);
        #endif
    }

    MPI_Bcast(types, num_particles, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(masses, num_particles, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(positions, num_particles, vectorMPI, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(velocities, chunk, vectorMPI, my_velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
}

void Output_State(double time, double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk) {
    int part;

    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
    if (my_rank == MASTER) {
        printf("Current time:%.2f\n", time);
        for (part = 0; part < num_particles; part++) {
            //       printf("%.3f ", masses[part]);
            printf("Particle:%3d\tX:%10.3e ", part, positions[part][X]);
            printf("\tY:%10.3e ", positions[part][Y]);
            printf("\tVx:%10.3e ", velocities[part][X]);
            printf("\tVy:%10.3e\n", velocities[part][Y]);
        }
        printf("\n");
    }
}

void update_state(unsigned char* image, int types[], double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk) {
    int part;

    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
    if (my_rank == MASTER) {
        for (part = 0; part < num_particles; part++) {
            update_img(image, positions[part][X], positions[part][Y], types[part], img_width, img_height);
        }
    }
}

void Compute_Force(int my_particles, double masses[], vector my_forces[], vector positions[], int num_particles, int chunk){

    int k, part;
    double m_g;
    vector f_part_k;
    double len, len_3, fact;

    /* Indice corrispondente alla particelle locali */
    part = my_rank*chunk + my_particles;
    my_forces[my_particles][X] = my_forces[my_particles][Y] = 0.0;

    #ifdef DEBUG_FORCES_BEFORE
    printf("Proc %d > Current total force on part %d = (%.3e, %.3e)\n", my_rank, part, my_forces[my_particles][X], my_forces[my_particles][Y]);
    #endif

    for (k = 0; k < num_particles; k++) {
        if (k != part) {
            if (masses[k] <= 0) {
                continue;
            }
            f_part_k[X] = positions[part][X] - positions[k][X];
            f_part_k[Y] = positions[part][Y] - positions[k][Y];


            len=sqrt(pow(f_part_k[X],2)+pow(f_part_k[Y],2));
            len_3=pow(len,3);


            m_g = G*masses[part]*masses[k];
            fact = m_g/len_3;

            f_part_k[X] *= fact;
            f_part_k[Y] *= fact;

            #ifdef DEBUG_FORCES_AFTER
            printf("Proc %d > Force on part %d due to part %d = (%.3e, %.3e)\n", my_rank, part, k, f_part_k[X], f_part_k[Y]);
            #endif

            /* Forza totale sulla particella */
            my_forces[my_particles][X] += f_part_k[X];
            my_forces[my_particles][Y] += f_part_k[Y];
        }
    }

}

void Update_Particles(int my_particles, double masses[], vector my_forces[], vector my_positions[], vector my_velocities[], int num_particles, int chunk, double delta_t) {

    int part;
    double fact;

    part = my_rank*chunk + my_particles;
    fact = delta_t/masses[part];

    #ifdef DEBUG_UPDATE_BEFORE
    printf("   Proc %d > Before update of %d:\n", my_rank, part);
    printf("   Position  = (%.3e, %.3e)\n", my_positions[my_particles][X], my_positions[my_particles][Y]);
    printf("   Velocity  = (%.3e, %.3e)\n", my_positions[my_particles][X], my_positions[my_particles][Y]);
    printf("   Net force = (%.3e, %.3e)\n", my_forces[my_particles][X], my_forces[my_particles][Y]);
    #endif

    my_positions[my_particles][X] += delta_t * my_velocities[my_particles][X];
    my_positions[my_particles][Y] += delta_t * my_velocities[my_particles][Y];
    my_velocities[my_particles][X] += fact * my_forces[my_particles][X];
    my_velocities[my_particles][Y] += fact * my_forces[my_particles][Y];

    #ifdef DEBUG_UPDATE_AFTER
    printf("Proc %d > Position of %d = (%.3e, %.3e), Velocity = (%.3e,%.3e)\n", my_rank, part, my_positions[my_particles][X], my_positions[my_particles][Y],
           my_velocities[my_particles][X], my_velocities[my_particles][Y]);
    #endif
}

void Generate_Output_File(double masses[], vector positions[], vector my_velocities[], int num_particles, int chunk){
    int line;
    MPI_Gather(my_velocities, chunk, vectorMPI, velocities, chunk, vectorMPI, MASTER, MPI_COMM_WORLD);
    if(my_rank==MASTER){
        FILE *fp_generate_output = fopen("generate_output.txt","w+");
        if(!fp_generate_output){
            printf("Cannot create output file <%s>\n", "generate_output");
            exit(1);
        }

        for(line=0; line<num_particles; line++){
            fprintf(fp_generate_output, "%lf %lf %lf %lf %lf\n", masses[line], positions[line][X], positions[line][Y], velocities[line][X], velocities[line][Y]);
        }

        fclose(fp_generate_output);

    }
}