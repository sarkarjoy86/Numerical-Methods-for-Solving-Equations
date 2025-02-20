/**
    Project Name: Numerical Methods
    ================================
    This program will calculate different Equations using differenet Numeric Methods.

    Prepared By:
    --------------------------------------------------
     Joy Sarkar                    - 0822220105101086
     Shawon Roy                    - 0822220105101075
     Mst. Farzana Akter Anamika    - 0822220205101094
    --------------------------------------------------
**/
// NB:The chrono Library is not properly working on The VS code as well as The Code Blocks, 
// it is prefered to run the code in Online compiler to View the Accurate Time Frame of each Methods. 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono> // For time calculations
#include <algorithm>
using namespace std;

// Function prototypes for Non Linear methods
double f(double x);
double f_prime(double x);
double g(double x);
double g_prime(double x);
double h(double x);
double h_prime(double x);
double i(double x);
double i_prime(double x);

// Functions for Jacobi and Gauss Seidal Method
#define f1_case1(x, y, z) (85 - 6 * y + z) / 27
#define f2_case1(x, y, z) (72 - 6 * x - 3 * z) / 15
#define f3_case1(x, y, z) (110 - x - y) / 54

#define f1_case2(x, y, z) (12.6 - 2 * y + 2 * z) / 26
#define f2_case2(x, y, z) (-14.3 - 3 * x - z) / 27
#define f3_case2(x, y, z) (6.01 - 2 * x - 3 * y) / 17

// Linear Method
void gaussElimination(double matrix[10][11], double result[10], int n);
void gaussJordan(double matrix[10][11], double result[10], int n);
void jacobiMethod(int choice);
void gaussSeidelMethod(int choice);

// Non - Linear Method
double bisection(double a, double b, int max_iterations, double (*func)(double));
double RegulaFalsi(double a, double b, int max_iterations, double (*func)(double));
double NewtonRaphson(double x, int max_iterations, double (*func)(double), double (*func_prime)(double));
double OnePointIteration(double x, int max_iterations, double (*func)(double));

int main()
{
    char problem_type;
    int method_choice, equationChoice;
    do
    {
        cout << endl
             << endl;
        cout << "\t\t\t\t --------------------------------------- " << endl;
        cout << "\t\t\t\t | Welcome to Numerical Methods System |" << endl;
        cout << "\t\t\t\t --------------------------------------- " << endl;

        cout << endl;
        cout << "\t\t\t\t\t    ---Main Menu---" << endl;
        cout << "\t\t\t\t    1. Solve Linear Equations" << endl;
        cout << "\t\t\t\t    2. Solve Non-Linear Equations" << endl;
        cout << "\t\t\t\t    3. Exit" << endl;
        cout << endl
             << "\t\t\t\t --------------------------------------- " << endl;
        cout << endl
             << endl;
        cout << "\t\t\t    Enter your choice (1-3): ";
        cin >> problem_type;

        // Cant escape untill pressing 1,2 or 3
        if (problem_type != '1' && problem_type != '2' && problem_type != '3')
        {
            cout << "\t\t Invalid choice. Please enter a valid option (1, 2, or 3)." << endl;
        }
    } while (problem_type != '1' && problem_type != '2' && problem_type != '3');

    switch (problem_type)
    {
    case '1':
        system("CLS");
        cout << endl
             << endl;
        cout << "\t\t\t\t ----------------------------------------------- " << endl;
        cout << "\t\t\t\t | Choose the Method to Solve Linear Equations |" << endl;
        cout << "\t\t\t\t ----------------------------------------------- " << endl;

        cout << endl;
        cout << "\t\t\t\t    1. Gauss Elimination" << endl;
        cout << "\t\t\t\t    2. Gauss-Jordan Elimination" << endl;
        cout << "\t\t\t\t    3. Jacobi Iteration Method" << endl;
        cout << "\t\t\t\t    4. Gauss-Seidel Iteration Method" << endl;
        cout << endl
             << "\t\t\t\t ----------------------------------------------- " << endl;
        cout << endl
             << endl;
        cout << "\t\t\t    Enter your choice (1-4): ";
        cin >> method_choice;
        cout << endl
             << endl;
        switch (method_choice)
        {
        case 1:
        {
            cout << endl
                 << endl;
            system("CLS");
            cout << "\t\t\t ------ Using Gauss Elimination Method ------" << endl
                 << endl;
            int n;
            double matrix[10][11], result[10];
            cout << "Enter the number of equations: ";
            cin >> n;
            cout << "Enter the coefficients and constants (row-wise):" << endl;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    cin >> matrix[i][j];
                }
            }
            gaussElimination(matrix, result, n);

            cout << "The solution is:" << endl;
            for (int i = 0; i < n; i++)
            {
                cout << "x" << i + 1 << " = " << result[i] << endl;
            }
            break;
        }

        case 2:
        {
            cout << endl
                 << endl;
            system("CLS");
            cout << "\t\t\t ------ Gauss-Jordan Elimination Method ------" << endl
                 << endl;
            int n;
            double matrix[10][11], result[10];
            cout << "Enter the number of equations: ";
            cin >> n;
            cout << "Enter the coefficients and constants (row-wise):" << endl;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    cin >> matrix[i][j];
                }
            }
            gaussJordan(matrix, result, n);

            cout << "The solution is:" << endl;
            for (int i = 0; i < n; i++)
            {
                cout << "x" << i + 1 << " = " << static_cast<int>(round(result[i])) << endl;
            }
            break;
        }

        case 3:
            cout << endl
                 << endl;
            system("CLS");
            cout << "\t\t\t ------ Jacobi Iteration Method ------" << endl
                 << endl;
            cout << "\t\t Choose the set of linear equations:\n";
            cout << "1. 27x + 6y - z = 85\n   6x + 15y + 3z = 72\n   x + y + 54z = 110\n";
            cout << "2. 26x + 2y + 2z = 12.6\n   3x + 27y + z = -14.3\n   2x + 3y + 17z = 6.01\n";
            cout << endl
                 << "\t Enter your choice (1 or 2): ";
            cin >> equationChoice;
            jacobiMethod(equationChoice);
            break;

        case 4:
            cout << endl
                 << endl;
            system("CLS");
            cout << "\t\t\t ------ Gauss-Seidel Iteration Method ------" << endl
                 << endl;
            cout << "\t\t Choose the set of linear equations:\n";
            cout << "1. 27x + 6y - z = 85\n   6x + 15y + 3z = 72\n   x + y + 54z = 110\n";
            cout << "2. 26x + 2y + 2z = 12.6\n   3x + 27y + z = -14.3\n   2x + 3y + 17z = 6.01\n";
            cout << endl
                 << "\t Enter your choice (1 or 2): ";
            cin >> equationChoice;
            gaussSeidelMethod(equationChoice);
            break;

        default:
            cout << "\t\t\t Invalid choice!! Please select a valid method." << endl;
            return 1;
        }
        break;

    case '2':
    {
        system("CLS");
        int equation_choice;
        cout << endl
             << endl;
        cout << "\t\t\t\t -------------------------------- " << endl;
        cout << "\t\t\t\t | Choose the Equation to Solve |" << endl;
        cout << "\t\t\t\t -------------------------------- " << endl;
        cout << endl;
        cout << "\t\t\t\t ----- Trigonometric Equations----" << endl;
        cout << "\t\t\t\t    1. 3x - cos(x) - 1" << endl;
        cout << "\t\t\t\t    2. 2*x + sin(x) - 3" << endl;
        cout << endl
             << "\t\t\t\t ----- Plynominal Equations----" << endl;
        cout << "\t\t\t\t    3. 2*x^3 + x^2 - 5*x - 1" << endl;
        cout << "\t\t\t\t    4. x^3 - 4*x - 9" << endl;
        cout << endl
             << "\t\t\t\t -------------------------------- " << endl;
        cout << endl
             << endl;
        cout << "\t\t\t Enter your choice (1-4): ";
        cin >> equation_choice;
        cout << endl
             << endl;
        double (*selected_function)(double);       // Functions
        double (*selected_function_prime)(double); // Functions for Derrivation
        double (*g_function)(double);              // g(x) for One point Iterations
        double initial_guess;

        switch (equation_choice)
        {
        case 1:
            selected_function = f;
            selected_function_prime = f_prime;
            // For One point Iteration
            initial_guess = 3.5;
            g_function = [](double x)
            { return (cos(x) + 1.0) / 3.0; }; // Defination of the function by Lambda Expression
            break;

        case 2:
            selected_function = g;
            selected_function_prime = g_prime;
            // For One point Iteration
            initial_guess = 3.5;
            g_function = [](double x)
            { return (3 - sin(x)) / 2.0; }; // Defination of the function by Lambda Expression
            break;

        case 3:
            selected_function = h;
            selected_function_prime = h_prime;
            // For One point Iteration
            initial_guess = 0.5;
            g_function = [](double x)
            { return cbrt((2 * x + 5) / 2.0); }; // Defination of the function by Lambda Expression
            break;

        case 4:
            selected_function = i;
            selected_function_prime = i_prime;
            // For One point Iteration
            initial_guess = 0.5;
            g_function = [](double x)
            { return cbrt(4 * x + 9); }; // Defination of the function by Lambda Expression
            break;

        default:
            cout << "\t\t\t Invalid choice!! Please select a valid equation (1-4)." << endl;
            return 1;
        }

        double a, b;
        int max_iterations = 100;

        cout << "Enter the Lower Guess (a): ";
        cin >> a;
        cout << "Enter the Upper Guess (b): ";
        cin >> b;

        // Measure runtime for Bisection Method
        auto start = chrono::high_resolution_clock::now(); // Calculates Starting Timne
        double root_bisection = bisection(a, b, max_iterations, selected_function);
        auto end = chrono::high_resolution_clock::now();                 // Calculates Ending Timne
        chrono::duration<double, micro> elapsed_bisection = end - start; // Calculates Timne Duration

        // Measure runtime for Regula Falsi Method
        start = chrono::high_resolution_clock::now();
        double root_regula_falsi = RegulaFalsi(a, b, max_iterations, selected_function);
        end = chrono::high_resolution_clock::now();
        chrono::duration<double, micro> elapsed_regula_falsi = end - start;

        // Measure runtime for Newton-Raphson Method
        start = chrono::high_resolution_clock::now();
        double root_newton_raphson = NewtonRaphson(a, max_iterations, selected_function, selected_function_prime);
        end = chrono::high_resolution_clock::now();
        chrono::duration<double, micro> elapsed_newton_raphson = end - start;

        // Measure runtime for One-Point Iteration Method
        start = chrono::high_resolution_clock::now();
        double root_one_point = OnePointIteration(initial_guess, max_iterations, g_function);
        end = chrono::high_resolution_clock::now();
        chrono::duration<double, micro> elapsed_one_point = end - start;

        // Printing results and runtimes
        cout << "Root found using Bisection method: " << fixed << setprecision(3) << root_bisection << endl;
        cout << "Time taken by Bisection method: " << fixed << setprecision(3) << elapsed_bisection.count() << " microseconds" << endl;
        cout << endl;
        cout << "Root found using Regula Falsi method: " << fixed << setprecision(3) << root_regula_falsi << endl;
        cout << "Time taken by Regula Falsi method: " << fixed << setprecision(3) << elapsed_regula_falsi.count() << " microseconds" << endl;
        cout << endl;
        cout << "Root found using Newton-Raphson method: " << fixed << setprecision(3) << root_newton_raphson << endl;
        cout << "Time taken by Newton-Raphson method: " << fixed << setprecision(3) << elapsed_newton_raphson.count() << " microseconds" << endl;
        cout << endl;
        cout << "Root found using One-Point Iteration method: " << fixed << setprecision(3) << root_one_point << endl;
        cout << "Time taken by One-Point Iteration method: " << fixed << setprecision(3) << elapsed_one_point.count() << " microseconds" << endl;
        cout << endl;

        // Find the method with the lowest runtime
        double min_time = min({elapsed_bisection.count(), elapsed_regula_falsi.count(),
                               elapsed_newton_raphson.count(), elapsed_one_point.count()});
        string best_method;
        // Checking the Best Method
        if (min_time == elapsed_bisection.count())
            best_method = "Bisection";
        else if (min_time == elapsed_regula_falsi.count())
            best_method = "Regula Falsi";
        else if (min_time == elapsed_newton_raphson.count())
            best_method = "Newton-Raphson";
        else if (min_time == elapsed_one_point.count())
            best_method = "One-Point Iteration";

        cout << "The fastest method is - " << best_method << " with a time of " << fixed << setprecision(3) << min_time << " microseconds." << endl;
        break;
    }
    case '3':
        cout << "\t\t -------Thank you for using the Numerical System!! GoodBye !!-------" << endl;
        break;
    default:
        cout << "Invalid problem type. Please select a valid option." << endl;
        return 1;
    }

    return 0;
}
// Define the functions and their derivatives
double f(double x)
{
    return (3.0 * x - cos(x) - 1.0);
}

double f_prime(double x)
{
    return (3.0 + sin(x));
}

double g(double x)
{
    return (2 * x + sin(x) - 3);
}

double g_prime(double x)
{
    return (2 + cos(x));
}

double h(double x)
{
    return (2 * pow(x, 3) - 2 * x - 5);
}

double h_prime(double x)
{
    return (6 * pow(x, 2) - 2);
}

double i(double x)
{
    return (pow(x, 3) - 4 * x - 9);
}

double i_prime(double x)
{
    return (3 * pow(x, 2) - 4);
}

// Bisection method implementation
double bisection(double a, double b, int max_iterations, double (*func)(double))
{
    if (func(a) * func(b) >= 0)
    {
        cout << "Bisection method cannot guarantee convergence within the given interval." << endl;
        return -1;
    }

    double c = 0;
    for (int i = 0; i < max_iterations; i++)
    {
        c = (a + b) / 2.00;
        if (func(c) == 0.0)
            break;
        if (func(c) * func(a) < 0)
            b = c;
        else
            a = c;
    }
    return c;
}

// Regula Falsi method implementation
double RegulaFalsi(double a, double b, int max_iterations, double (*func)(double))
{
    if (func(a) * func(b) >= 0)
    {
        cout << "Regula-Falsi method cannot guarantee convergence within the given interval." << endl;
        return -1;
    }

    double c = 0;
    for (int i = 0; i < max_iterations; i++)
    {
        c = ((a * func(b)) - (b * func(a))) / (func(b) - func(a));
        if (func(c) == 0.0)
            break;
        if (func(c) * func(a) < 0)
            b = c;
        else if (func(c) * func(b) < 0)
            a = c;
    }
    return c;
}

// Newton-Raphson method implementation
double NewtonRaphson(double x, int max_iterations, double (*func)(double), double (*func_prime)(double))
{
    for (int i = 0; i < max_iterations; ++i)
    {
        double fx = func(x);
        double fpx = func_prime(x);
        if (fabs(fpx) < 1e-8)
        {
            cout << "Derivative near zero; Newton-Raphson method may fail." << endl;
            return x;
        }

        double x_new = x - fx / fpx;
        if (fabs(x_new - x) < 0.01)
            return x_new;

        x = x_new;
    }

    cout << "Newton-Raphson method did not converge within the maximum number of iterations." << endl;
    return x;
}

// One-Point Iteration method implementation
double OnePointIteration(double initial_guess, int max_iterations, double (*func)(double))
{
    double x = initial_guess;
    double x_prev = x; // Variable to store previous iterate
    for (int i = 0; i < max_iterations; i++)
    {
        double gx = func(x); // Compute g(x)
        x = gx;
        // To stop the iteration
        if (fabs(x - x_prev) < 0.001) // Check if difference between consecutive x values is small
        {
            return x;
        }
        x_prev = x; // Update x_prev for next iteration
    }
    return x;
}

// Linear Eguations
// // Function to perform Gauss elimination
void gaussElimination(double matrix[10][11], double result[10], int n)
{
    int stepCount = 0; // Initialize step counter

    // Forward elimination process
    for (int i = 0; i < n; i++)
    {
        // Partial pivoting: find the maximum element in the current column
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i]))
            {
                maxRow = k;
            }
        }

        // Swap the maximum row with the current row
        for (int k = 0; k < n + 1; k++)
        {
            swap(matrix[maxRow][k], matrix[i][k]);
        }
        stepCount++; // Count the swap as a step

        // Make the diagonal element 1 (normalization)
        double diag = matrix[i][i];
        for (int j = i; j < n + 1; j++)
        {
            matrix[i][j] /= diag;
        }
        stepCount++; // Count the normalization as a step

        // Eliminate the elements below the pivot
        for (int k = i + 1; k < n; k++)
        {
            double factor = matrix[k][i];
            for (int j = i; j < n + 1; j++)
            {
                matrix[k][j] -= factor * matrix[i][j];
            }
            stepCount++; // Count each row elimination as a step
        }
    }

    // Back substitution process
    for (int i = n - 1; i >= 0; i--)
    {
        result[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++)
        {
            result[i] -= matrix[i][j] * result[j];
        }
        stepCount++; // Count each step in the back substitution as a step
    }

    // Output the total step count
    cout << endl
         << "Total steps needed: " << stepCount << endl;
}

// Function to perform Gauss-Jordan elimination
void gaussJordan(double matrix[10][11], double result[10], int n)
{
    int stepCount = 0; // Initialize step counter

    // Forward elimination process
    for (int i = 0; i < n; i++)
    {
        // Partial pivoting: find the maximum element in the current column
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(matrix[k][i]) > fabs(matrix[maxRow][i]))
            {
                maxRow = k;
            }
        }

        // Swap the maximum row with the current row
        for (int k = 0; k < n + 1; k++)
        {
            swap(matrix[maxRow][k], matrix[i][k]);
        }
        stepCount++; // Count the swap as a step

        // Normalize the current row
        double diag = matrix[i][i];
        if (fabs(diag) < 1e-10)
        {
            cout << "Matrix is singular or nearly singular." << endl;
            return;
        }
        for (int j = 0; j < n + 1; j++)
        {
            matrix[i][j] /= diag;
        }
        stepCount++; // Count the normalization as a step

        // Eliminate all other rows
        for (int k = 0; k < n; k++)
        {
            if (k != i)
            {
                double factor = matrix[k][i];
                for (int j = 0; j < n + 1; j++)
                {
                    matrix[k][j] -= factor * matrix[i][j];
                }
                stepCount++; // Count the row elimination as a step
            }
        }
    }

    // Output the result and the step count
    cout << endl
         << "Total steps needed: " << stepCount << endl;

    // Extract the results from the matrix
    for (int i = 0; i < n; i++)
    {
        result[i] = matrix[i][n];
    }
}

// Function to perform Jacobi Iteration Method
void jacobiMethod(int choice)
{
    const float EPSILON = 0.001;
    float x0 = 0, y0 = 0, z0 = 0, x1, y1, z1, e1, e2, e3;
    int step = 1;
    bool continueIteration = true;

    cout << fixed << setprecision(3);

    do
    {
        switch (choice)
        {
        case 1:
            x1 = f1_case1(x0, y0, z0);
            y1 = f2_case1(x0, y0, z0);
            z1 = f3_case1(x0, y0, z0);
            break;
        case 2:
            x1 = f1_case2(x0, y0, z0);
            y1 = f2_case2(x0, y0, z0);
            z1 = f3_case2(x0, y0, z0);
            break;
        default:
            cout << "Invalid equation choice!" << endl;
            return;
        }

        e1 = fabs(x0 - x1);
        e2 = fabs(y0 - y1);
        e3 = fabs(z0 - z1);
        // To check 2 values of x,y and z are same or not
        if (e1 <= EPSILON && e2 <= EPSILON && e3 <= EPSILON)
        {
            continueIteration = false;
        }
        // Updating
        x0 = x1;
        y0 = y1;
        z0 = z1;
        step++;

    } while (continueIteration);

    cout << endl
         << "Solution: x = " << x1 << ", y = " << y1 << " and z = " << z1 << endl;
    cout << endl
         << "Total steps needed to converge: " << step - 1 << endl;
    cout << endl
         << "\t -- Try to Use Gauss - Seidal Method to reduce the step Number -- " << step - 1 << endl;
    cout << "\t --And more Steps means a more Time Complexity -- " << endl;
}

// Function to perform Gauss Seidal Method
void gaussSeidelMethod(int choice)
{
    const float EPSILON = 0.001; // Define the epsilon value
    float x0 = 0, y0 = 0, z0 = 0, x1, y1, z1, e1, e2, e3;
    int step = 1;
    bool continueIteration = true;

    cout << fixed << setprecision(3);

    do
    {
        switch (choice)
        {
        case 1:
            x1 = f1_case1(x0, y0, z0);
            y1 = f2_case1(x1, y0, z0); // Use the updated x1 for y1 calculation
            z1 = f3_case1(x1, y1, z0); // Use the updated x1 and y1 for z1 calculation
            break;
        case 2:
            x1 = f1_case2(x0, y0, z0);
            y1 = f2_case2(x1, y0, z0); // Use the updated x1 for y1 calculation
            z1 = f3_case2(x1, y1, z0); // Use the updated x1 and y1 for z1 calculation
            break;
        default:
            cout << "Invalid equation choice!" << endl;
            return;
        }

        e1 = fabs(x0 - x1);
        e2 = fabs(y0 - y1);
        e3 = fabs(z0 - z1);

        if (e1 <= EPSILON && e2 <= EPSILON && e3 <= EPSILON)
        {
            continueIteration = false;
        }
        x0 = x1;
        y0 = y1;
        z0 = z1;
        step++;

    } while (continueIteration);
    cout << endl
         << "Solution: x = " << x1 << ", y = " << y1 << " and z = " << z1 << endl;
    cout << endl
         << "Total steps needed to converge: " << step - 1 << endl;
    cout << endl
         << "\t --Gauss - Seidal Iteration Method Needs fewer steps than Jacobi Iteration Method-- " << endl;
    cout << "\t --And Fewer Steps means a less Time Complexity -- " << endl;
}