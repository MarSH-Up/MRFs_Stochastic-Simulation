/* Assigment 3: MRFs Optimization
 * Description: Implement the stochastic simulation algorithm for obtaining the most probable configuration of
 *              a 1st order regular MRF considering the at least two variantes ICM and Metropolis, using MAP
 *
 * Version: 1.0
 * Date: 2021 March
 * By: Mario De Los Santos
 * Instituto Nacional de Astrofísica Óptica y Electónica.
 *
 * Class: Probabilistic Graphical Models
 * By: PhD Enrique Sucar
 *
 * Reference; Sucar, L. E. (2015). Probabilistic graphical models.
 *            Advances in Computer Vision and Pattern Recognition.
 *            London: Springer London. doi, 10(978), 2
 */
#include <iostream>
#include<math.h>
using namespace std;
/*Version Notes:
 * Right now the project should return the matrix, both we are displaying it directly in the functions as a void
 * - Still working in two versions:
 *                       - Version 2.0 would include the matrix return to the main function, so you can operate with it.
 *                       - Version 3.0, would include the use of the project with binary images.
 *                       Date; March 3rd
 * */

/*Class Description:
 * Version 1.0
 *
 * The class is composed of:
 *                          - 1 Parametrizated constructor
 *                          - 2 return functions (float, int)
 *                          - 2 Operational fuctions
 *                          - 1 Destructor (Empty in the version 1)
 *  Each one is used in the class, is not recommended to delete any of them without confirming the application in the other functions.
 *
 *  Base operation of the class:
 *      1. Object creation:
 *          You would need to create an object with some parameters
 *              Example:  MRF_Optimization T1(px, py, (int*)Observation_matrix);
 *      2. Public fuctions:
 *          The fuctions to call and use the MRF operations need some parameters also. Those depends on the user exercise:
 *              Example: Following the example in the point 1.
 *                     -  T1.Stochastic_ICM_MAP(lamda,No_interactions);
 *                     -  T1.Stochastic_ICM_Metropolis(lamda,No_interactions,Metropolis_Prob);
 *
 * Inside the class you would find the documentation for each part, and some relevant information about the code.
 * */
class MRF_Optimization
{
    /*Variables used along the code, the idea is to have a predefined size,
     * if you are going to work with images I would recommend you to delimete the size with the maximum resolution that you would accept
     * Example:
     *
     *     #define default_size_x 1920
     *     #define default_size_y 1080
     *
     * Which is the resolution of pixels in a fullHD image.
     */
    #define default_size_x 10
    #define default_size_y 10
    //Delimitation variables for the incoming matrixes
    int px, py;
    //Lamda value needed to calculate the energies of the "pixeles"
    //Reference: Chapter 6, pag 99. 6.4 Inference. Book in the head of the code
    int lamda_clss;
    //Number of iterations that you want to run the method selected, you can determinate any number,
    //the program would let you know when the matrix stop changing
    int iterations;

    //The following matrices are to keep the class private, adherents to the data hiding philosophy
    int Observation_mtx[default_size_x][default_size_y]={0};
    int MRF_F[default_size_x][default_size_y]={0};
    int MRF_F1[default_size_x][default_size_y]={0};

public:
    /*The constructor would need 3 parameters:
     *  (Size in X of the matrix, size in Y, and the observation matrix(the binary image matrix if you are using images))
     *  The call is described as an example in the class head
     */
    MRF_Optimization(int size_px, int size_py, int * Obs)
    {
        //We confirm that the sizes are diferente of 0, so we can continue
        if(size_px&&size_py !=0)
        {
            //Data hiding
            px=size_px;
            py=size_py;
            //Data hiding, adding the observation matrix to the object information
            for(int i=0;i<px;i++){
                for(int j=0;j<py;j++) {
                    cout<<(Observation_mtx[i][j] = *(Obs + i * py + j));
                }cout<<endl;}
            //Just a warranty to get the matrixes clean with the 1 and 0 respectively
            for(int i=0;i<default_size_x;i++)
                for(int j=0;j<default_size_y;j++) {
                    MRF_F1[i][j]=1;
                    MRF_F[i][j]=0;
                }
        }
        //If the parameters are wrong, your would be out
        else
            cout<<"Error with the parametrization"<<endl;
        //Deleting the pointer
        delete Obs;
    }

    /* The call is described below, we request the lamda value to calculate the energies,
     * and the number of iterections desired
     *
     *  Example.Stochastic_IMC_MAP(4, 1); //Example
     * */
    void Stochastic_ICM_MAP(int lamda, int intr)
    {
        cout<<"ICM_MAP"<<endl;
        iterations  = intr;
        lamda_clss = lamda;
        float z0=0,z1=0;
        int iteraction_review=0;
        int aux[default_size_x][default_size_y]={0};
        //Define a aux matrix to work with in the interactions
        for(int i=0;i<px;i++)
            for (int j = 0; j < py; j++)
                aux[i][j] = MRF_F[i][j];


        /*  "There are several variants of this general algorithm according to variations on different aspects. One is the way in which the optimal configuration is defined, for
        * which there are two main alternatives: MAP and MPM. In the Maximum A posteriori Probability or MAP, the optimum configuration is taken as the configuration at the
        * end of the iterative process. In the case of Maximum Posterior Marginals or MPM, the most frequent value for each variable in all the iterations is taken as the optimum
        * configuration." -Book reference
        *
        * Also, consider that the book provides a better pseudocode with the approach, we use it as reference here.
        *  Loop that would run the iteractions called by the user
        *  */
        //In MAP usually would take 1 iteraction to reconstruct the observation matrix
        for(int s=0;s<iterations +1;s++)
        {
            //Run around the rows
            for(int i=0;i<px;i++)
            {
                z0 = 0; //Rest Z after each iteraction, the function of Z0 and Z1 would be explained in the  Smoothing functions
                z1 = 0;
                //Run across the columns
                for(int j=0;j<py;j++)
                {
                    //We calculate the energies for each case the 0 and 1 case for the MRF matrix in comparison with the Observation one
                    z0 = Smoothing((int*)aux,i,j);
                    z1 = Smoothing((int*)MRF_F1,i,j);
                    cout<<z0<<" "<<z1<<endl; //Just to confirm the decision

                    //Base selection, we would take the option with less energy.
                    /* Example: We take the following observation Matrix:
                     *
                     * Observation = { 0 0 0 0
                     *                 0 1 1 0
                     *                 0 1 0 0
                     *                 0 0 1 0}
                     *
                     * Consider that our MRF matrix would be a 4x4 matrix with ceros, we calculate the energies as explained in the smoothing fuction and compare
                     *
                     *  Consider MFR[1][1] as 0 and for the Observation Matrix as 1, then the energies would be:
                     *          0 = 4 and 1 = 3, so considering that 1 has less energy, we remplace that position in the MRF original matrix.
                     *
                     *   The same step would be repited the number of iteractions desired, or until the matrix doesnt chance anymore.
                     *
                     * */
                    if(z0>=z1)
                        MRF_F[i][j] = MRF_F1[i][j];
                }
            }
            //The matrix check function, would basically confirm if the two matrix are the same. If is the case, we would have a 1 in the variable
            iteraction_review = Matrix_Check(MRF_F,Observation_mtx);

            // If Iteraction review is 1, we return the number of iteractions needed ant then brake the iteractions loop
            if(iteraction_review == 1) {
                cout<<"Finished in: "<<s<<" interactions"<<endl;
                break;}
        }
      //For this version we use the vizualization just to confirm the results, the future development in this code is review in the head of the code
      cout<<"Vizualization, just debbuging"<<endl;
        for(int i=0;i<px;i++) {
            for (int j = 0; j < py; j++) {
                cout<<(MRF_F[i][j]);}
        cout<<" "<<endl;
        }
    }
    void Stochastic_ICM_Metropolis(int lamda, int intr, float P)
    {
        cout<<"Metropolis"<<endl;
        iterations  = intr;
        lamda_clss = lamda;
        float z0=0,z1=0;
        int iteraction_review=0;
        int aux[default_size_x][default_size_y]={0};

        //Define a aux matrix to work with in the interactions
        for(int i=0;i<px;i++)
            for (int j = 0; j < py; j++)
                aux[i][j] = MRF_F[i][j];

        //Run around the iteractions
        for(int s=0;s<iterations ;s++)
        {
            //Run across the rows
            for(int i=0;i<px;i++)
            {
                z0 = 0; //Rest Z after each iteraction, the function of Z0 and Z1 would be explained in the  Smoothing functions
                z1 = 0;
                //Run across the columns
                for(int j=0;j<py;j++)
                {
                    //We calculate the energies for each case the 0 and 1 case for the MRF matrix in comparison with the Observation one
                    z0 = Smoothing((int*)aux,i,j);
                    z1 = Smoothing((int*)MRF_F1,i,j);                    //Base selection, we would take the option with less energy.
                    /* Example: We take the following observation Matrix:
                     *
                     * Observation = { 0 0 0 0
                     *                 0 1 1 0
                     *                 0 1 0 0
                     *                 0 0 1 0}
                     *
                     * Consider that our MRF matrix would be a 4x4 matrix with ceros, we calculate the energies as explained in the smoothing fuction and compare
                     *
                     *  Consider MFR[1][1] as 0 and for the Observation Matrix as 1, then the energies would be:
                     *          0 = 4 and 1 = 3 So, considering that in the case of Metropolis we consider two things, first that the Z0 energy is bigger or equal to the Z1,
                     *                          The we generate a random number, between 0 and 1, the idea is that if the random number pass the probaility given by the user,
                     *                          then we change the 0 to 1, if not, we keep it as 0.
                     *
                     *   The same step would be repited the number of iteractions desired, or until the matrix doesnt chance anymore.
                     *
                     * */

                    if(z0>=z1) {
                        float uf = float(rand() % (11027 + 1)) / 11027;
                        cout<<uf<<endl;
                        if(uf>P)
                            MRF_F[i][j]=MRF_F1[i][j];
                        /*
                         * "Metropolis: with a fixed probability, P, it selects a configuration with a higher energy. "
                         * "ICM provides the most efficient approach, and can be used if the optimization landscape is known to be convex or relatively “simple”. SA is the most robust approach
                         *  against local minima, nevertheless it tends to be the most costly computationall. Metropolis provides a middle ground between ICM and SA. The choice of which
                         *  variation to use depends on the application, in particular the size of the model complexity of the landscape and efficiency considerations..."
                         *
                         * */
                    }
                }
            }
            //The matrix check function, would basically confirm if the two matrix are the same. If is the case, we would have a 1 in the variable
            iteraction_review = Matrix_Check(MRF_F,Observation_mtx);

            // If Iteraction review is 1, we return the number of iteractions needed ant then brake the iteractions loop
            if(iteraction_review == 1) {
                cout<<"Finished in: "<<s<<" interactions"<<endl;
                break;}
        }
        //For this version we use the vizualization just to confirm the results, the future development in this code is review in the head of the code
        for(int i=0;i<px;i++) {
            for (int j = 0; j < py; j++) {
                cout << MRF_F[i][j];
                MRF_F[i][j] = 0;
            }
            cout<<endl;
        }
    }
    ~MRF_Optimization(){}

private:
    float Smoothing(int *MRF, int row, int colum)
    {
        /* Energy Fuction:
         *  Considering a regular MRF of order n, the energy functio can be expressed in terms of functions of subsets of completely connected variables of different sizes, 1, 2, 3, ...:
         *  Given the Gibbs equivalence, the problem of finding the configuration of maximum probability for a MRF is transformed to finding the configuration of minimum energy
         *  Then the following code follow the next fuction:
         *      Up(f) = EcVc(f) + Lambda * EoVo(f)
         *
         *   So we obtain the values of each pixel (consider a pixel each point in the matrix) making a compromise between both constaints:
         *
         *      Vc(f) = (f - u)^2
         *      Vo(f) = (f - g)^2
         *
         *    Where g is the observation matrix, and f de MRF matrix.
         *
         *    Example: Consider the matrix F = 0 0 0 0  and G = 0 0 0 0 and a Lambda of 4
         *                                     0 0 0 0          0 1 1 0
         *                                     0 0 0 0          0 1 0 0
         *                                     0 0 0 0          0 0 1 0
         *      Then we are going to operate the point (0,0) in F:
         *          So, Vo = (0-0)^2 + (0 - 0)^2 + Lambda * (0-0)^2
         *          The easier way to understand it is:
         *              Vo = (MRF[i][j] - MRF[i+1][j])^2 + (MRF[i][j] - MRF[i][j+1])^2 + Lambda * (MRF[i][j]*Observation[i][j])^2
         * */
        int MRF_I[default_size_x][default_size_y]={0};
        for(int i=0;i < px;i++) {
            for (int j = 0; j < py; j++)
                MRF_I[i][j] = *(MRF + i * py + j);
        }
        float z=0;
        if (row == 0) {
            if (colum == 0) {
                //Upper left case
                z = abs(pow(((MRF_I[row][colum]) - (MRF_F[1][colum])), 2) +
                        pow(((MRF_I[row][colum]) - (MRF_F[row][colum + 1])), 2) +
                        (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
            }
            else if (colum == (py - 1)) //Upper Right case
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row+1][colum]), 2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum-1]), 2)  + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
            else //Any case between the two already mentioned
                z = abs( pow((MRF_I[row][colum])- (MRF_F[row+1][colum]),2) + pow((MRF_I[row][colum]) - (MRF_I[row][colum+1]),2) +
                         pow((MRF_I[row][colum]) - (MRF_F[row][colum-1]),2) + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
        }
        else if(row == (px-1)){
            if(colum==0) {//Lower left case
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row - 1][colum]), 2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum + 1]), 2) +
                        (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
                /*cout << (MRF_I[row][colum]) - (MRF_F[row - 1][colum]) << "+"
                     << (MRF_I[row][colum]) - (MRF_F[row][colum + 1]) << "+"
                     << (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)) << endl;*/
            }
            else if (colum == (py - 1)) //Lower right case
                z = abs(pow((MRF_I[row][colum])-MRF_F[row-1][colum],2) + pow((MRF_I[row][colum])-(MRF_F[row][colum-1]),2) +
                        (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
            else //Any case between the two already mentioned
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row+1][colum]),2) + pow((MRF_I[row][colum]) - (MRF_F[row][colum+1]),2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum-1]),2) + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
        }
        else { //All the cases in the matrix that are not a coner, basically
            if(colum==0) //If the point is in the left edge
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row + 1][colum]), 2) + pow((MRF_I[row][colum]) - (MRF_F[row - 1][colum]), 2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum+1]), 2) + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
            else if(colum == (py-1)) //If the market is in the right  edge
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row + 1][colum]), 2) + pow((MRF_I[row][colum]) - (MRF_F[row - 1][colum]), 2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum - 1]), 2) + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
            else //All the cases between, basically a criss-cross
                z = abs(pow((MRF_I[row][colum]) - (MRF_F[row + 1][colum]), 2) + pow((MRF_I[row][colum]) - (MRF_F[row - 1][colum]), 2) +
                        pow((MRF_I[row][colum]) - (MRF_F[row][colum + 1]), 2) + pow((MRF_I[row][colum]) - (MRF_F[row][colum - 1]), 2) + (lamda_clss * pow((MRF_I[row][colum]) - (Observation_mtx[row][colum]), 2)));
        }
        return z;
        delete MRF;
    }
    int Matrix_Check(int A[][default_size_y],int B[][default_size_y])
    {
        //Basic comparation between matrixes, here a link if you want to check it:
        // https://www.geeksforgeeks.org/program-to-check-if-two-given-matrices-are-identical/
        for(int i=0;i<px;i++)
            for(int j=0;j<py;j++)
                if(A[i][j]!=B[i][j])
                    return 0;

        return 1;
    }

};
int main() {
    //The following staments are just examples of how to declare the requirement to operate the basic model of the class
    //This using matrixes as example and testing.
    #define px 4 //Matrix's sizee
    #define py 4//Matrix's size

    int Observation_matrix[px][py]={{0,0,0,0}, //Observation Matrix for testing
                                    {0,1,1,0},
                                    {0,1,0,0},
                                    {0,0,1,0}};
    int lamda = 4; //Lamba constant needed to calculate the energies
    int No_interactions=10; //Number of interactions for both process ICM MAP and ICM Metropolis
    float Metropolis_Prob = 0.99997; //Probability value

    MRF_Optimization T1(px, py, (int*)Observation_matrix); //Object creation, read the description in the class

    /*Fuctions calls, both needs the lamda and intereactions to operate
     * To more detail on how to use them review the class documentation
     */
    T1.Stochastic_ICM_MAP(lamda,No_interactions);
    //T1.Stochastic_ICM_Metropolis(lamda,No_interactions,Metropolis_Prob);
    return 0;
}
