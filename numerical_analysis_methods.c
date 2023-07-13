#include <stdio.h>
#include <math.h>

# define MAX 100

// Structures
typedef struct Polynom {
    int degree;
    double coefficients[MAX];
}poly;

// Helping Functions
void mainmenu();
poly readPolynom(poly);
int readMatrix(double [MAX][MAX]);
void printMatrix(double [][MAX], int);
void printMatrixGE(double [][MAX], int);
void printArray(double [], int);
poly derivative(poly);
double polynomOutput(poly, double);

// Main Functions
double bisectionMethod();
double regulafalsiMethod();
double newtonraphsonMethod();
int inverseMatrix();
int gaussElimination();
int gaussSeidal();
double numericalDerivative();
double trapezMethod();
double simpsonMethod();
double gregoryNewtonInterpolation();


void main() {
    mainmenu();
}

void mainmenu(){
    int choice = 11;

    while(choice != 0){
        printf("\n---------------------------------\n");
        printf("Quit: 0\nBisection: 1\nRegula-Falsi: 2\nNewton-Raphson: 3\nInverse of a NxN matrix: 4\nGauss Elemination: 5\nGauss Seidal: 6\nNumerical Derivative: 7\nSimpson: 8\nTrapez: 9\nGregory Newton Interpolation: 10\nYour Choice: ");
        scanf("%d", &choice);
        printf("---------------------------------\n");

        if(choice == 1) bisectionMethod();
        if(choice == 2) regulafalsiMethod();
        if(choice == 3) newtonraphsonMethod();
        if(choice == 4) inverseMatrix();
        if(choice == 5) gaussElimination();
        if(choice == 6) gaussSeidal();
        if(choice == 7) numericalDerivative();
        if(choice == 8) simpsonMethod();
        if(choice == 9) trapezMethod();
        if(choice == 10) gregoryNewtonInterpolation();
    }
}

// Read Polynom
poly readPolynom(poly p){
    int i;
    
    printf("Enter the degree of the polynom: ");
    scanf("%d", &p.degree);
    for ( i = 0; i < p.degree+1; i++){
        printf("Enter coefficient of x^%d: ", p.degree-i);
        scanf("%lf", &p.coefficients[p.degree-i]);
    }

    // fonksiyonu manuel olarak girme   
    //p.degree = 4; p.coefficients[4] = 1; p.coefficients[3] = 0; p.coefficients[2] = -3; p.coefficients[1] = 7; p.coefficients[0] = 1;

    // Print Function
    printf("your function is: ");
    for ( i = 0; i < p.degree+1; i++){
        if(i != p.degree){
            if(p.coefficients[p.degree-i] == 0);
            else if(p.coefficients[p.degree-i] == 1)printf("x^%d + ", p.degree-i);
            else printf("%lfx^%d + ", p.coefficients[p.degree-i], p.degree-i);
        }
        else printf("%lf", p.coefficients[p.degree-i]);
    }
    return p;
}

// Result of an x value for a polynom
double polynomOutput(poly p, double x){
    int i;
    double result = 0;
    for ( i = 0; i < p.degree+1; i++) result += p.coefficients[i]*pow(x, i);
    return result;
}

// Analytical Polynom Derivative
poly derivative(poly p){
    int i;
    for(i=0; i < p.degree; i++){
        p.coefficients[i] = p.coefficients[i+1]*(i+1);
    }
    p.degree--;
    return p;
}

// Read Matrix (it also returns the size of the matrix)
int readMatrix(double m[MAX][MAX]){
    int i,j,n;
    /*double deneme[10][10] ={{3,4,-3},
                             {3,1,-2},
                             {1,-1,4}};*/
    printf("\nEnter the size of the matrix: ");
    scanf("%d", &n);  // n = 3; 

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            printf("[%d,%d]: ", i, j);
            scanf("%lf", &m[i][j]); //m[i][j] = deneme[i][j];
        }
    }
    printf("\nYour matrix is: ");
    printMatrix(m, n);
    return n;
}

void printMatrix(double m[][MAX], int n){
    int i, j;
    for (i = 0; i < n; i++){
        printf("\n");
        for (j = 0; j < n; j++){
            printf("%lf ", m[i][j]);
        }
    }
}

void printMatrixGE(double m[][MAX], int n){
    int i, j;
    for (i = 0; i < n; i++){
        printf("\n");
        for (j = 0; j < n+1; j++){
            if (j == n) printf("= ");
            printf("%lf ", m[i][j]);
        }
    }
}

void printArray(double arr[], int n){
    int i;
    for (i = 0; i < n; i++){
        printf("x%d = %lf ", i+1, arr[i]);
    }
}

// 1) Bisection Method
double bisectionMethod(){
    printf("\nYou chose Bisection Method. \n");
    poly p;
    double  start, end, mid;
    double epsilon, error;
    int i = 0, stopcriteria = 0, quit = 0, maxiteration ;

    p = readPolynom(p);

    printf("\nChoose a Stopping Criteria:\n    f(x) <= epsilon: 1\n    absolute Error <= epsilon: 2\n");
    scanf("%d", &stopcriteria); //stopcriteria = 1;
    if((stopcriteria != 1) && (stopcriteria != 2)){
        printf("\nPlease choose an existing criteria. ");
        return 1;
    }

    printf("\nEnter the range(ex:5 7): ");
    scanf("%lf %lf", &start, &end); //start = 2; end = 3;    
    // mean value theorem
    if( polynomOutput(p, start)*polynomOutput(p, end) > 0 ){
        printf("\nPlease choose another range. ");
        return 1;
    }

    printf("\nEnter an epsilon value: ");
    scanf("%lf", &epsilon); //epsilon = .01;

    printf("\nEnter Maximum iteration number: ");
    scanf("%d", &maxiteration); //maxiteration = 50;

    while((quit != 1) && (i < maxiteration)){
        mid = (start+end)/2.;
        if (polynomOutput(p, mid)*polynomOutput(p, end) < .0) start = mid;
        else end = mid;
        if(stopcriteria == 1){
            error = fabs(polynomOutput(p, mid));
            if(error <= epsilon) quit = 1;
        }
        else {
            error = (end-start)/pow(2., i);
            if( error <= epsilon) quit = 1;
        }
        printf("\n%d. iteration: f(%lf) = %lf\nError = %lf",i+1,mid, polynomOutput(p, mid), error);
        i++;
    }
    printf("\n===================================");
    printf("\nRESULT: %d. iteration: f(%lf) = %lf\n",i,mid, polynomOutput(p, mid));
    printf("\n===================================");
    return mid;
}

// 2) Regula-Falsi Method
double regulafalsiMethod(){
    printf("\nYou chose Regula-Falsi Method. \n");
    poly p;
    double  start, end, mid;
    double epsilon, error;
    int i = 0, stopcriteria, quit = 0, maxiteration;

    p = readPolynom(p);

    printf("\nChoose a Stopping Criteria:\n    f(x)<=epsilon: 1\n    absolute Error<=epsilon: 2\n");
    scanf("%d", &stopcriteria);
    if((stopcriteria != 1) && (stopcriteria != 2)){
        printf("\nPlease choose an existing criteria. ");
        return 1;
    }

    printf("\nEnter the range(ex:5 7): ");
    scanf("%lf %lf", &start, &end); //start = 2; end = 3;

    if( polynomOutput(p, start)*polynomOutput(p, end) > 0 ){
        printf("\nPlease choose another range. ");
        return 1;
    }

    printf("\nEnter an epsilon value: ");
    scanf("%lf", &epsilon); // epsilon = .01;

    printf("\nEnter Maximum iteration number: ");
    scanf("%d", &maxiteration); //maxiteration = 50;

    while((quit != 1) && (i < maxiteration)){
        mid = (end*polynomOutput(p, start)-start*polynomOutput(p, end))/(polynomOutput(p, start)-polynomOutput(p, end));
        if (polynomOutput(p, mid)*polynomOutput(p, end) < .0) start = mid;
        else end = mid;
        if(stopcriteria == 1){
            error = fabs(polynomOutput(p, mid));
            if(error <= epsilon) quit = 1;
        }
        else {
            error = (end-start)/pow(2., i);
            if( error <= epsilon) quit = 1;
        }
        printf("\n%d. iteration: f(%lf) = %lf\nError = %lf",i+1,mid, polynomOutput(p, mid), error);
        i++;
    }
    printf("\n===================================");
    printf("\nRESULT: %d. iteration: f(%lf) = %lf",i,mid, polynomOutput(p, mid));
    printf("\n===================================");
    return mid;
}

// 3) Newton-Raphson Method
double newtonraphsonMethod(){
    printf("\nYou chose Newton Raphson Method. \n");
    poly p;
    double  x1, x0, difference, der = 1;
    double epsilon;
    int i = 0, maxiteration;

    p = readPolynom(p);

    printf("\nEnter an x value to start: ");
    scanf("%lf", &x0); //x0 = 1;

    printf("\nEnter an epsilon value: ");
    scanf("%lf", &epsilon); //epsilon = 0.001;
    difference = epsilon+1;

    printf("\nEnter Maximum iteration number: ");
    scanf("%d", &maxiteration); maxiteration = 50; 

    while((difference > epsilon) && (i < maxiteration) && (der != 0)){
        der = polynomOutput(derivative(p), x0);
        x1 = x0 - (polynomOutput(p, x0)/der);
        difference = fabs(x1-x0);
        printf("\n%d. iteration: f(%lf) = %lf\ndifference: %lf",i+1,x0, polynomOutput(p, x0), difference);
        x0 = x1;
        i++;
    }
    if(der != 0){
        printf("\n===================================");
        printf("\n%d. iteration: RESULT = %lf",i, x0);
        printf("\n===================================");
    }
    else printf("\nNewton-Raphson Method is not suitable for this function.");
    return x0;
}

// 4) Inverse of a NxN matris
int inverseMatrix(){
    printf("\nYou chose Inverse Matrix. \n");
    int i, j, k, n, zero;
    double tmp1, tmp2;
    double m[MAX][MAX];
    double inverse[MAX][MAX] = {{1,0,0},{0,1,0},{0,0,1}};

    n = readMatrix(m);

    for(k = 0; k < n; k++){
        // Controlling horizontally
        zero = 0;
        while ((m[k][zero] == 0)) zero++;

        // If main one of the diagonal elements is 0 change the line
        if (zero < n) zero = 0;
        while (((m[k][k] == 0) && (zero < n))){
            // Vertical
            zero++;
            for(i = 0; i < n; i++){
                if((k+zero) < n){
                    tmp2 = m[k][i];
                    m[k][i] = m[k+zero][i];
                    m[k+zero][i] = tmp2;

                    tmp2 = inverse[k][i];
                    inverse[k][i] = inverse[k+zero][i];
                    inverse[k+zero][i] = tmp2;
                }
                else{
                    tmp2 = m[k][i];
                    m[k][i] = m[0][i];
                    m[0][i] = tmp2;

                    tmp2 = inverse[k][i];
                    inverse[k][i] = inverse[0][i];
                    inverse[0][i] = tmp2;
                } 
            }
        }

        if(zero < n){
            tmp1 = m[k][k];
            for (i = 0; i < n; i++){
                // Divide the line by the main diagonal element of the line
                inverse[k][i] = inverse[k][i]/tmp1;
                m[k][i] = m[k][i]/tmp1;
            }
            for (i = 0; i < n; i++){
                tmp2 = m[i][k];
                for (j = 0; j < n; j++){
                    if (i != k){
                        inverse[i][j] += -tmp2*inverse[k][j];
                        m[i][j] += -tmp2*m[k][j];
                    }
                }
            }
            
            // Print the matrix
            printf("\n\nStep %d: ", k+1);
            printf("\n=======================");
            printf("\nFirst Matrix: ");
            printMatrix(m, n);
            printf("\n\nInverse matrix: ");
            printMatrix(inverse, n);
        }
        else{
            printf("\nThis Matris is not Inversible\n");
            return 1;
        }
        printf("\n=======================");
    }

}

// 5) Gauss Elimination of a NxN matris
int gaussElimination(){
    printf("\nYou chose Gauss Elimination Method. \n");
    int i, j, k, n, zero;
    double tmp1, tmp2;
    double m[MAX][MAX];
    double results[MAX];
    //double arr[MAX] = {-8,9,1}; // Bu silinecek

    printf("\nPlease enter coefficients of the variables as a matrix: ");
    n = readMatrix(m);

    printf("\n\nPlease enter the values that each row equals: ");
    for ( i = 0; i < n; i++) scanf("%lf", &m[i][n]); //for ( i = 0; i < n; i++)  m[i][n] = arr[i];

    for(k = 0; k < n; k++){
        // Controlling horizontally
        zero = 0;
        while ((m[k][zero] == 0)) zero++;

        // If one of the main diagonal elements is 0, change the line
        if (zero < n) zero = 0;
        while (((m[k][k] == 0) && (zero < n))){
            // Vertical
            zero++;
            for(i = 0; i < n;i++){
                if((k+zero) < n){
                    tmp2 = m[k][i];
                    m[k][i] = m[k+zero][i];
                    m[k+zero][i] = tmp2;
                }
                else{
                    tmp2 = m[k][i];
                    m[k][i] = m[0][i];
                    m[0][i] = tmp2;
                } 
            }
        }
        if(zero < n){
            tmp1 = m[k][k];
            m[k][n] = m[k][n]/tmp1;
            for (i = 0; i < n; i++){
                m[k][i] = m[k][i]/tmp1;
            }
            for (i = k; i < n; i++){
                tmp1 = m[i][k];
                for (j = 0; j < n+1; j++){
                    if (i != k){
                        m[i][j] += -tmp1*m[k][j];
                    }
                }
            }
            // Print the steps
            printf("\n\nStep %d: ", k+1);
            printf("\n=======================\n");
            printMatrixGE(m, n);
            printf("\n=======================\n");
            
            // Finding values by going backwards
            for (i = n-1; i >= 0; i--){
                results[i] = 0;
                for ( j = i; j < n; j++) results[i] += results[j]*m[i][j];
                results[i] = m[i][n] - results[i];
            }
        }
        else{
            printf("Wrong Matris. try again\n");
            return 1;
        }
    }
    if(k == n){
        printf("\n=======================\nResults: ");
        for ( i = 0; i < n; i++) printf("x%d = %lf  ", i+1, results[i]);
        printf("\n=======================\n");
    }
}

// 6) Gauss Seidal Iteration
int gaussSeidal(){ 
    printf("\nYou chose Gauss-Seidal Method. \n");
    int i, j, k, n, zero, max, maxiteration, end = 0;
    double tmp1, x, epsilon;
    double m[MAX][MAX];
    double results[MAX] = {0};
    //double deneme[MAX] = {-8,9,1}; 

    printf("\nPlease enter coefficients of the variables as a matrix: ");
    n = readMatrix(m);

    // Controlling horizontally and vertically
    for ( i = 0; i < n; i++){
        zero = 0;
        while ((m[i][zero] == 0)) zero++;
        if(zero < n){
            zero = 0;
            while ((m[i][zero] == 0)) zero++;
        }
    }
    
    if(zero < n){
        printf("\n\nPlease enter the values that each row equals: ");
        for ( i = 0; i < n; i++) scanf("%lf", &m[i][n]); //for ( i = 0; i < n; i++)  m[i][n] = deneme[i];

        printf("\nEnter the maximum iteration number: ");
        scanf("%d", &maxiteration); //maxiteration = 50;

        printf("\nEnter the epsilon value: ");
        scanf("%lf", &epsilon); //epsilon = .001;

        // Put the highest values to main diagonal
        for (i = 0; i < n; i++){
            max = i;
            for (j = i+1; j < n; j++){
                if(fabs(m[j][i]) > fabs(m[max][i])) max = j;
            }
            for(j = 0; j < n+1; j++){
                tmp1 = m[i][j];
                m[i][j] = m[max][j];
                m[max][j] = tmp1;
            }
        }

        k = 0;
        while(!end && (k < maxiteration)){
            end = 1;
            for (i = 0; i < n; i++){
                x = 0;
                for (j = 0; j < n; j++){
                    if (j != i) x += -m[i][j]*results[j];
                    else x += m[i][n];
                }
                // controlling the error
                if(fabs(results[i]-(x/m[i][i])) > epsilon) end = 0;
                results[i] = x/m[i][i];
            }
            printf("\nIteration %d: \n", k+1);
            printArray(results, n);
            k++;
        }
        if (k == maxiteration)
            printf("\nThe results diverged. Gauss-Seidal is not working for this matrix. ");
    }
    else {
        printf("Wrong Matrix. Try again");
        return 1;
    }
}
    
// 7) Numerical Derivative
double numericalDerivative(){
    printf("\nYou chose Numerical Derivative. \n");
    poly p;
    double  x, result;
    double delta;
    int choice;

    p = readPolynom(p);

    printf("\nChoose one of them: \n    Forward: 1 \n    Backward: 2 \n    Central: 3\n");
    scanf("%d", &choice); //choice = 3;

    printf("\nEnter a x value: ");
    scanf("%lf", &x); //x = 1;

    printf("\nEnter a delta value: ");
    scanf("%lf", &delta); //delta = .1;

    if(choice == 1) result = ((polynomOutput(p, x) - polynomOutput(p, x-delta))/delta);
    else if(choice == 2) result = ((polynomOutput(p, x+delta) - polynomOutput(p, x))/delta);
    else result = ((polynomOutput(p, x+delta) - polynomOutput(p, x-delta))/(2.*delta));
    
    printf("\n===================================");
    printf("\n*Analytical Derivative of f(%lf) = %lf \nNumerical Derivative of f(%lf) = %lf",x, polynomOutput(derivative(p), x),x, result);
    printf("\n===================================");
    return result;
}

// 8) Simpson Method
double simpsonMethod(){
    printf("\nYou chose Simpson Methods. \n");
    poly p;
    double start, end, h, result = 0;
    int i, n, choice;

    p = readPolynom(p);

    printf("\nChoose one of them: \n    1/3 method: 1 \n    3/8 method: 2\n");
    scanf("%d", &choice); //choice = 2;

    printf("\nEnter a range(ex: 1 2): ");
    scanf("%lf %lf", &start, &end); //start = 0; end = 6;

    printf("\nEnter a n value: ");
    scanf("%d", &n); //n = 1;

    

    // Simpson 1/3
    if (choice == 1){
        h = (end-start)/(2*n);
        result = polynomOutput(p, start) + polynomOutput(p, end);
        for (i = 1; i < n*2; i++){   
            if(i%2==0) result += 2.*polynomOutput(p, start+h*i);
            else result += 4.*polynomOutput(p, start+h*i);
        }
        result = h*result/3;
    }
    // Simpson 3/8
    if (choice == 2){
        h = (end-start)/(3*n);
        result = polynomOutput(p, start) + polynomOutput(p, end) ;
        for (i = 1; i < 3*n; i++){   
            if(i%3==0) result += 2*polynomOutput(p, start+h*i);
            else result += 3*polynomOutput(p, start+h*i);
        }
        result = h*result*(3./8.);
    }
    
    printf("\n===================================");
    printf("\nIntegral [%lf to %lf] = %lf by Simpson 1/3 Method",start, end, result);
    printf("\n===================================");
    return result;
}

// 9) Trapez Method
double trapezMethod(){
    printf("\nYou chose Trapez Method. \n");
    poly p;
    double start, end, h, result = 0;
    int i, n;

    p = readPolynom(p);

    printf("\nEnter a range(ex: 1 2): ");
    scanf("%lf %lf", &start, &end); //start = 0; end = 4; 

    printf("\nEnter a n value: ");
    scanf("%d", &n); //n = 5;

    h = (end-start)/n;

    for (i = 0; i < n; i++){
        result += h/2*(polynomOutput(p, start+h*i)+polynomOutput(p, start+h*(i+1)));
    }
    
    printf("\n===================================");
    printf("\nIntegral [%lf to %lf] = %lf by Trapez Method",start, end, result);
    printf("\n===================================");
    return result;
}

// 10) Gregory Newton Enterpolasyonu
double gregoryNewtonInterpolation(){ 
    printf("\nYou chose Gregory-Newton Method. \n");
    int i, j, n, degree = 0;
    double k, x, h, x0, answer;
    double results[MAX];
    double deltafx[MAX];

    printf("\nEnter the first x value: ");
    scanf("%lf", &x0); //x0 = 2;

    printf("\nEnter the increase in x values: ");
    scanf("%lf", &h); // h = 2;

    printf("\nEnter the number of points: ");
    scanf("%d", &n); //n = 5;

    printf("\nEnter which value you want to estimate: ");
    scanf("%lf", &x); //x = 2;

    for (i = 0; i < n; i++){
        printf("\nEnter the value of f(%lf) = ", x0+i*h);
        scanf("%lf", &results[i]);
    }
    
    while (results[0] != 0){
        deltafx[degree] = results[0];
        for (i = 0; i < n-1; i++) results[i] = results[i+1] - results[i];
        n--;
        degree++;    
    }

    k = (x-x0)/h;
    answer = deltafx[0];
    for (i = 1; i < degree; i++){
        x0 = 1; // I am using this as a temp
        for (j = 0; j < i; j++)
           x0 *= (double)(k-j)/(double)(j+1.);
        answer += x0*deltafx[i];
    }
    printf("\n===================================");
    printf("\nThe Value you want: f(%lf) = %lf", x, answer);
    printf("\n===================================");
}