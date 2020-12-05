#include <iostream>
#include <iomanip>
#include <fstream>
#include <ios>
#include <list>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iterator>
#include<unordered_map>
#include<map>

using namespace std;

using namespace std::chrono;

int n = 0;
class Point {

public:
    double xval;
    double yval;
    Point()
    {

    }
    Point(const Point &other)
    {
        this->xval = other.xval;
        this->yval = other.yval;
    }
};


class GridPoint {
    int xval, yval;  // First and last names

    GridPoint(int x, int y)
    {
        xval = x;
        yval = y;
    }

    // Match both first and last names in case
    // of collisions.
    bool operator==(const GridPoint& p) const
    {
        return xval == p.xval && yval == p.yval;
    }
};


class MyHashFunction {
public:
    unsigned long long first, second;
    template <class T1, class T2>

    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
//        return hash1 ^ hash2;
          return hash1 + hash2 ;
    }
//    template <class T3, class T4>
//    bool operator==(const pair<T3, T4>& p) const
//    {
//        return first == p.first && second == p.second;
//    }
};


Point A;
Point B;
vector<Point> coords;
vector<Point> coordspart3;
vector<Point> coordspart4;
Point finalPoints2[2];
double part2distance;
double timePart2;
Point finalPoints1[2];
double part1distance;
double timePart1;
double min2 = 1;
double min3 = 1;
double timePart3;
Point finalPoints3[2];
double part3distance;

Point finalPoints4[2];
double timePart4;
double part4distance;


void part1();
void part2();
void part3();
double bruteforce(list<Point> list, int n); // part 1

double min(double x, double y);
double dist(Point point, Point point1);
void generatePoints();
void readPoints(vector<Point> &abc);
int compareX(const void *point, const void *point1);
bool compareYY(Point x, Point y);
double closestRecursive(vector<Point> &vector);
double closestRecursivepart3(vector<Point> &vectorPoints);
double stripPart2(vector<Point> &in, double d);
double stripPart3(vector<Point> &in, double d);
double bruteForce(vector<Point> &vector, int n); // part 2
double bruteForcePart3(vector<Point> &vector, int n); //part3
int minInt(int i, int i1);
list<Point> total;

void part4();

//map<Point,Point> umap;
std::unordered_map<pair<unsigned long long, unsigned long long>, Point, MyHashFunction> umap;
void randomize();



int main() {

    std::ofstream output;
    output.open ("results.txt");


    srand( time(NULL) );


//    part1();
//
//    cout << "Points with Brute Force Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints1[0].xval << " , " << finalPoints1[0].yval << ") , ("<< finalPoints1[1].xval << " , " << finalPoints1[1].yval << ")"
//         <<"\n"<<"Minimum Distance with Brute Force Approach: " << part1distance << "\n";
//
//    cout << "Time taken by Part 1: "
//         << timePart1 << " microseconds" << endl << "\n";

  //  part2();
//    cout << "N is: " << n << "\n";

//    cout << "Points with Recursive Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints2[0].xval << " , " << finalPoints2[0].yval << ") , ("<< finalPoints2[1].xval << " , " << finalPoints2[1].yval << ")"
//         <<"\n"<<"Minimum Distance with Recursive Approach: " << part2distance << "\n";
//
//    cout << "Time taken by Part 2: "
//         << timePart2 << " microseconds" << endl << "\n";
//
//    output << "Points with Recursive Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints2[0].xval << " , " << finalPoints2[0].yval << ") , ("<< finalPoints2[1].xval << " , " << finalPoints2[1].yval << ")"
//           <<"\n"<<"Minimum Distance with Recursive Approach: " << part2distance << "\n";
//
//    output << "Time taken by Part 2: "
//           << timePart2 << " microseconds" << endl << "\n";
    part3();



    output << "Points with Complete Recursive Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints3[0].xval << " , " << finalPoints3[0].yval << ") , ("<< finalPoints3[1].xval << " , " << finalPoints3[1].yval << ")"
           <<"\n"<<"Minimum Distance with Complete Recursive Approach: " << part3distance << "\n" ;

    output << "Time taken by Part 3: "
           << timePart3 << " microseconds" << endl << "\n";



    cout << "Points with Complete Recursive Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints3[0].xval << " , " << finalPoints3[0].yval << ") , ("<< finalPoints3[1].xval << " , " << finalPoints3[1].yval << ")"
           <<"\n"<<"Minimum Distance with Complete Recursive Approach: " << part3distance << "\n";

    cout << "Time taken by Part 3: "
           << timePart3 << " microseconds" << endl << "\n";

    part4();

    output << "Points with Randomized Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints4[0].xval << " , " << finalPoints4[0].yval << ") , ("<< finalPoints4[1].xval << " , " << finalPoints4[1].yval << ")"
           <<"\n"<<"Minimum Distance with Complete Recursive Approach: " << part4distance << "\n" ;

    output << "Time taken by Part 4: "
           << timePart4 << " microseconds" << endl;



    cout << "Points with Randomized Approach: " << fixed<< setprecision(23)<< "(" <<finalPoints4[0].xval << " , " << finalPoints4[0].yval << ") , ("<< finalPoints4[1].xval << " , " << finalPoints4[1].yval << ")"
           <<"\n"<<"Minimum Distance with Complete Recursive Approach: " << part4distance << "\n" ;

    cout << "Time taken by Part 4: "
           << timePart4 << " microseconds" << endl;

    cout << (part4distance == part3distance) << endl;

    output.close();
    return 0;
}

void part4() {
    readPoints(coordspart4);
//    for(Point p: coordspart4)
//    {
//        cout << "B4RanX: " << p.xval << "  B4RanY: " << p.yval << "\n";
//    }
    auto start = high_resolution_clock::now();
    randomize();
    bool firstTime = false;
    double theta = dist(coordspart4[0], coordspart4[1]);

   // cout << "theta: " << theta << endl;
    double d = theta/2;
  //  cout << "d: " << d << endl;
    //change how finding grid points
    unsigned long long x = coordspart4[0].xval/d;
    unsigned long long y = coordspart4[0].yval/d;
  //  cout << coordspart4[0].xval << "  " << coordspart4[0].yval << endl;

  //  cout << x << "  " << y << endl;
    //GridPoint p(x, y);
    pair<unsigned long long, unsigned long long> p(x,y);
    umap[p] = coordspart4[0];

    x = coordspart4[1].xval/d;
    y = coordspart4[1].yval/d;

    finalPoints4[0] = coordspart4[0];
    finalPoints4[1] = coordspart4[1];
    //GridPoint p2(x, y);
    pair<unsigned long long, unsigned long long> p2(x,y);
    umap[p2] = coordspart4[1];

    bool foundSmaller = false;

    for(int i = 2; i < (int)coordspart4.size(); i++)
    {

//        cout << "D: " << d << endl;

        //change how finding grid points
        unsigned long long h = coordspart4[i].xval/d;
        unsigned long long j = coordspart4[i].yval/d;

//        cout << "X: " << h << "Y: " << j << "\n";
//        cout << "Xact: " << coordspart4[i].xval << "Yact: " << coordspart4[i].yval << "\n";
        //GridPoint gen(x,y);
        pair<unsigned long long, unsigned long long> gen(h,j);
        unsigned long long leftSide = h-2;
        unsigned long long rightSide = h+2;
        unsigned long long topSide = j+2;
        unsigned long long bottomSide = j-2;

        if(leftSide < 0)
        {
            leftSide = 0;
        }
        if(bottomSide < 0)
        {
            bottomSide = 0;
        }

        if(rightSide > (unsigned long long)(1/d))
        {
            rightSide = (unsigned long long)(1/d);
        }
        if(topSide > (unsigned long long)(1/d))
        {
            topSide = (unsigned long long)(1/d);
        }
            // add if any are negative make 0
            // if any are greater than 1/d than make 1/d
        for(unsigned long long a = leftSide; a <= rightSide; a++ ) {
            for(unsigned long long b = bottomSide; b <= topSide; b++)
            {

                pair<unsigned long long, unsigned long long> check(a,b);
                if(a == rightSide)
                {
                    umap[gen] = coordspart4[i];
                }
                if(umap.find(check) != umap.end()) //!(check.first == gen.first && check.second == gen.second))
                {

                    if(((umap.at(check).xval != coordspart4[i].xval || umap.at(check).yval != coordspart4[i].yval))) {
//                        cout << "Here bud!" << endl;
//                        cout << check.first << "  " << check.second << endl;

                        if (dist(umap.at(check), coordspart4[i]) <= theta) {
                            theta = dist(umap.at(check), coordspart4[i]);
                            d = theta / 2;
                            foundSmaller = true;
                            finalPoints4[0] = umap.at(check);
                            finalPoints4[1] = coordspart4[i];

//                            cout << "distance less" << "  " << d << endl;
//
//                            cout << coordspart4[i].xval << "  " << coordspart4[i].yval << endl;
//                            cout << umap.at(check).xval << "  " << umap.at(check).yval << endl;

                        }


                           // umap[gen] = coordspart4[i];
                    }


                }

                firstTime = true;

            }


        }
        if(foundSmaller == true) {
            foundSmaller = false;
            //cout << "i: " << i << endl;
            umap.clear();
  //          cout << "CLEARED" << endl;
            for (int n = 0; n <= i; n++) {
                unsigned long long c = coordspart4[n].xval / d;
                unsigned long long e = coordspart4[n].yval / d;
//                cout << coordspart4[n].xval << "  " << coordspart4[n].yval << endl;
//                cout << c << "  " << e << endl;
                //GridPoint p(c, e);
                pair<unsigned long long, unsigned long long> p3(c, e);
                umap[p3] = coordspart4[n];
            }
          //  cout << "Size: " << umap.size() << endl;


        }


    }

   // cout << "distance: " << theta << endl;

    part4distance = theta;
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);


    timePart4 = (double)duration.count();

//    for(Point p: coordspart4)
//    {
//        cout << p.xval<< "  "<< p.yval << "\n";
//    }

}

void randomize() {
    srand (time(NULL));

    for (int i = n - 1; i > 0; i--)
    {
        int j = rand() % (i + 1);

        Point temp = coordspart4[i];
        coordspart4[i] = coordspart4[j];
        coordspart4[j] = temp;
    }
}

void part3()
{
    readPoints(coordspart3);
    auto start = high_resolution_clock::now();
    qsort(&coordspart3[0], n, sizeof(Point), compareX);
    part3distance = closestRecursivepart3(coordspart3);//, 0, coordspart3.size());
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);


    timePart3 = (double)duration.count();

}
double closestRecursivepart3(vector<Point> &vectorPoints) {

    if ((int)vectorPoints.size() <= 3) {
        return bruteForcePart3(vectorPoints, (int)vectorPoints.size());

    }

    int mid = (int)vectorPoints.size()/2;
    Point midPoint = vectorPoints[mid];

    std::vector<Point> split_lo(vectorPoints.begin(), vectorPoints.begin() + mid);
    std::vector<Point> split_hi(vectorPoints.begin() + mid + 1, vectorPoints.end());

    double dl = closestRecursivepart3(split_lo);
    double dr = closestRecursivepart3(split_hi);


    double d = min(dl, dr);



    vector<Point> strip;
    for (int i = 0; i < (int)vectorPoints.size(); i++) {
        if (abs(vectorPoints[i].xval - midPoint.xval) <= d)
            strip.push_back(vectorPoints[i]);
    }


    d = stripPart3(strip, d);


    return d;

}
double stripPart3(vector<Point> &strip, double d)
{

    sort(strip.begin(), strip.end(),  compareYY);
    for (int i = 0; i < (int)strip.size()-1; i++)
        for (int j = i + 1; j<minInt(i+15,(int)strip.size()); j++)
            if (dist (strip[i], strip[j]) < min3) {
                min3 = dist (strip[i], strip[j]);
                d = dist(strip[i], strip[j]);


                Point *part3Pointer;
                part3Pointer = new Point(strip[i]);
                finalPoints3[0] = *part3Pointer;

                part3Pointer = new Point(strip[j]);
                finalPoints3[1]= *part3Pointer;


            }
    return d;

}

int minInt(int i, int i1) {
    if(i<i1)
        return i;
    else return i1;
}

double bruteForcePart3(vector<Point> &vector, int n) {


    std::vector<Point>::iterator outside;
    for ( outside = vector.begin(); outside != vector.end(); outside++) {
        for (auto inside = next(outside, 1); inside != vector.end(); inside++){

            if ( dist(*outside, *inside) < min3) {

                min3 = dist(*outside, *inside);
                A = *outside;
                B = *inside;

                Point *part3Pointer;
                part3Pointer = new Point(A);
                finalPoints3[0] = *part3Pointer;

                part3Pointer = new Point(B);
                finalPoints3[1]= *part3Pointer;

            }
        }
    }


    double d = dist(A, B);
    return d;
}
void generatePoints() {
    std::ofstream myfile;
    myfile.open ("points.txt");

    for(int i = 0; i < 500; i++)
    {
        Point generate;
        generate.xval =  ((double) rand() / (RAND_MAX));
        generate.yval =  ((double) rand() / (RAND_MAX));

        total.push_back(generate);
        myfile << std::fixed << std::setprecision(23) << generate.xval << "  " << generate.yval<< endl;

    }

    myfile.close();
}

void part2() {
    readPoints(coords);


    auto start = high_resolution_clock::now();




    qsort(&coords[0], n, sizeof(Point), compareX);


    part2distance = closestRecursive(coords);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);

    timePart2 = duration.count();

}

double closestRecursive(vector<Point> &vectorPoints) {

    if ((int)vectorPoints.size() <= 3) {

        return bruteForce(vectorPoints, (int)vectorPoints.size());
    }

    int mid = (int)vectorPoints.size()/2;
    Point midPoint = vectorPoints[mid];

    std::vector<Point> split_lo(vectorPoints.begin(), vectorPoints.begin() + mid);
    std::vector<Point> split_hi(vectorPoints.begin() + mid + 1, vectorPoints.end());


    double dl = closestRecursive(split_lo);
    double dr = closestRecursive(split_hi);


    double d = min(dl, dr);


    vector<Point> strip;
    for (int i = 0; i < (int)vectorPoints.size(); i++) {
        if (abs(vectorPoints[i].xval - midPoint.xval) <= d)
            strip.push_back(vectorPoints[i]);
    }

    d = stripPart2(strip, d);

    return d;

}

double bruteForce(vector<Point> &vector, int n) {


    std::vector<Point>::iterator outside;
    for (outside = vector.begin(); outside != vector.end(); outside++) {
        for (auto inside = next(outside, 1); inside != vector.end(); inside++){
           // double dis = dist(*outside, *inside);
            if ( dist(*outside, *inside) < min2) {
                min2 = dist(*outside, *inside);
                A = *outside;
                B = *inside;

                Point *part2Pointer;
                part2Pointer = new Point(A);
                finalPoints2[0] = *part2Pointer;

                part2Pointer = new Point(B);
                finalPoints2[1]= *part2Pointer;
            }
        }
    }

    double d = dist(A, B);




    return d;
}

double stripPart2(vector<Point> &strip, double d) {


    for (int i = 0; i < (int)strip.size(); ++i)
        for (int j = i + 1; j < (int)strip.size() ; j++)
            if (dist (strip[i], strip[j]) < min2) {
                min2 = dist (strip[i], strip[j]);
                d = dist(strip[i], strip[j]);

                Point *part2Pointer;
                part2Pointer = new Point(strip[i]);
                finalPoints2[0] = *part2Pointer;

                part2Pointer = new Point(strip[j]);
                finalPoints2[1]= *part2Pointer;
            }

    return d;
}


double min(double x, double y)
{
    if(x < y)
        return x;
    else
    return y;
}

void readPoints(vector<Point> &input) {


    ifstream inputFile("points.txt");

    double x,y;
    int w =0;

    while (inputFile>>x>>y) {

        Point generate;
        generate.xval = x;
        generate.yval = y;
        input.push_back(generate);
        w++;

    }
    if(n == 0)
    {
        n = w;
    }

inputFile.close();

}

int compareX(const void *a, const void *b) {
    const Point* x = (Point *) a;
    const Point* y = (Point *) b;

    if (x->xval > y->xval)
        return 1;
    else if (x->xval < y->xval)
        return -1;

    return 0;
}

bool compareYY(Point x, Point y) {

    if (x.yval > y.yval)
        return true;
    else if (x.yval < y.yval)
        return false;

    return 0;
}

void part1() {
   //cout << "in part 1" << "\n";

    generatePoints();

    auto start = high_resolution_clock::now();
    part1distance = bruteforce(total, n);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);

   // cout << "Time taken by Part 1: " << duration.count() << " microseconds" << endl;
    // cout << "Smallest Distance Part 1: " << fixed << setprecision(23) <<part1distance << "\n";

    timePart1 = duration.count();

}

double bruteforce(list<Point> list, int n) {

    double min = 1;
    std::list<Point>::iterator outside;
    for (outside = list.begin(); outside != list.end(); ++outside) {
        for (auto inside = next(outside, 1); inside != list.end(); ++inside){
            double dis = dist(*outside, *inside);
           if ( dis < min) {
                min = dis;
                A = *outside;
                B = *inside;
            }
        }
    }

    double d = dist(A, B);

    finalPoints1[0] = A;
    finalPoints1[1] =B;

    return d;
}

double dist(Point p1, Point p2) {
    return sqrt( pow((p1.xval - p2.xval), 2) +
                         pow((p1.yval - p2.yval), 2));
}


