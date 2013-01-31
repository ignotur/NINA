#include <fstream>
#include "stars.h"

using namespace std;

struct decl {
    double quant;
};

TMap::TMap (void) {

    size = 1038961;
    ifstream in_Tb  ("TbGal_tot_CasA1pix.bin", ios::binary);
    decl Tb_tmp;

    for (int i = 0; i < size; i++)	{

        in_Tb.seekg(sizeof(Tb_tmp)*i, ios_base::beg);
        in_Tb.read((char*) &Tb_tmp, sizeof(Tb_tmp));
        Tb[i] = Tb_tmp.quant;
    }
}

double TMap::get_Tb (int N) {
    return Tb[N];
}
