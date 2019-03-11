#include <iostream>

using namespace std;

extern "C" {
double d_to_dm (double, double, double); 
}

int main () {

double gl, gb, dist, res;

gl = 0.0;
gb = 0.0;
dist = 500.0;

cout << d_to_dm (gl, gb, dist) << endl;

//cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;
cout << gl << "\t" << gb << "\t" << dist << endl;

res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);
res= d_to_dm (gl, gb, dist);

cout << res << endl;

return 0;
}
