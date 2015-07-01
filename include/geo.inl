#include "global.hpp"
#include "geo.hpp"

vector<point> geo::getLocSpts(int eType, int order)
{
  vector<point> outPts;

  if (eType == TRI) {
    outPts.resize(2*order+1);
    if (order == 0) {
      outPts[0].x =  0.;
      outPts[0].y =  0.;
    }
    else if (order == 1) {
      outPts[0].x = -1.;
      outPts[1].x =  1.;
      outPts[2].x = -1.;

      outPts[0].y = -1.;
      outPts[1].y = -1.;
      outPts[2].y =  1.;
    }
    else if (order == 2) {
      outPts[0].x = -1.;
      outPts[1].x =  0.;
      outPts[2].x =  1.;
      outPts[3].x = -1.;
      outPts[4].x =  0.;
      outPts[5].x = -1.;

      outPts[0].y = -1.;
      outPts[1].y = -1.;
      outPts[2].y = -1.;
      outPts[3].y =  0.;
      outPts[4].y =  0.;
      outPts[5].y =  1.;
    }
    else if (order == 3) {
      outPts[0].x = -1.;
      outPts[1].x = -0.447213595499958;
      outPts[2].x =  0.447213595499958;
      outPts[3].x =  1.;
      outPts[4].x = -1.;
      outPts[5].x = -0.333333333333333;
      outPts[6].x =  0.447213595499958;
      outPts[7].x = -1.000000000000000;
      outPts[8].x = -0.447213595499958;
      outPts[9].x = -1.000000000000000;

      outPts[0].y = -1.;
      outPts[1].y = -1.;
      outPts[2].y = -1.;
      outPts[3].y = -1.;
      outPts[4].y = -0.447213595499958;
      outPts[5].y = -0.333333333333333;
      outPts[6].y = -0.447213595499958;
      outPts[7].y =  0.447213595499958;
      outPts[8].y =  0.447213595499958;
      outPts[9].y =  1.;
    }
    else if (order == 4) {
      outPts[0].x  = -1.000000000000000;
      outPts[1].x  = -0.654653670707977;
      outPts[2].x  =  0.000000000000000;
      outPts[3].x  =  0.654653670707977;
      outPts[4].x  =  1.000000000000000;
      outPts[5].x  = -1.000000000000000;
      outPts[6].x  = -0.551583507555306;
      outPts[7].x  =  0.103167015110611;
      outPts[8].x  =  0.654653670707977;
      outPts[9].x  = -1.000000000000000;
      outPts[10].x = -0.551583507555306;
      outPts[11].x =  0.000000000000000;
      outPts[12].x = -1.000000000000000;
      outPts[13].x = -0.654653670707977;
      outPts[14].x = -1.000000000000000;

      outPts[0].y  = -1.000000000000000;
      outPts[1].y  = -1.000000000000000;
      outPts[2].y  = -1.000000000000000;
      outPts[3].y  = -1.000000000000000;
      outPts[4].y  = -1.000000000000000;
      outPts[5].y  = -0.654653670707977;
      outPts[6].y  = -0.551583507555306;
      outPts[7].y  = -0.551583507555306;
      outPts[8].y  = -0.654653670707977;
      outPts[9].y  =  0.000000000000000;
      outPts[10].y =  0.103167015110611;
      outPts[11].y =  0.000000000000000;
      outPts[12].y =  0.654653670707977;
      outPts[13].y =  0.654653670707977;
      outPts[14].y =  1.000000000000000;
    }
    else if (order == 5) {
      outPts[0].x  = -1.000000000000000;
      outPts[1].x  = -0.765055323929465;
      outPts[2].x  = -0.285231516480645;
      outPts[3].x  =  0.285231516480645;
      outPts[4].x  =  0.765055323929465;
      outPts[5].x  =  1.000000000000000;
      outPts[6].x  = -1.000000000000000;
      outPts[7].x  = -0.684472514501909;
      outPts[8].x  = -0.171245477332075;
      outPts[9].x  =  0.368945029003817;
      outPts[10].x =  0.765055323929465;
      outPts[11].x = -1.000000000000000;
      outPts[12].x = -0.657509045335851;
      outPts[13].x = -0.171245477332074;
      outPts[14].x =  0.285231516480645;
      outPts[15].x = -1.000000000000000;
      outPts[16].x = -0.684472514501908;
      outPts[17].x = -0.285231516480645;
      outPts[18].x = -1.000000000000000;
      outPts[19].x = -0.765055323929465;
      outPts[20].x = -1.000000000000000;

      outPts[0].y  = -1.000000000000000;
      outPts[1].y  = -1.000000000000000;
      outPts[2].y  = -1.000000000000000;
      outPts[3].y  = -1.000000000000000;
      outPts[4].y  = -1.000000000000000;
      outPts[5].y  = -1.000000000000000;
      outPts[6].y  = -0.765055323929465;
      outPts[7].y  = -0.684472514501909;
      outPts[8].y  = -0.657509045335851;
      outPts[9].y  = -0.684472514501909;
      outPts[10].y = -0.765055323929465;
      outPts[11].y = -0.285231516480645;
      outPts[12].y = -0.171245477332074;
      outPts[13].y = -0.171245477332074;
      outPts[14].y = -0.285231516480645;
      outPts[15].y =  0.285231516480645;
      outPts[16].y =  0.368945029003817;
      outPts[17].y =  0.285231516480645;
      outPts[18].y =  0.765055323929465;
      outPts[19].y =  0.765055323929465;
      outPts[20].y =  1.000000000000000;
    }
    else if (order == 6) {
      outPts[0].x =  -1.000000000000000;
      outPts[1].x =  -0.830223896278567;
      outPts[2].x =  -0.468848793470714;
      outPts[3].x =   0.000000000000000;
      outPts[4].x =   0.468848793470714;
      outPts[5].x =   0.830223896278567;
      outPts[6].x =   1.000000000000000;
      outPts[7].x =  -1.000000000000000;
      outPts[8].x =  -0.769516496193442;
      outPts[9].x =  -0.364769256536689;
      outPts[10].x =  0.108864652571010;
      outPts[11].x =  0.539032992386883;
      outPts[12].x =  0.830223896278567;
      outPts[13].x = -1.000000000000000;
      outPts[14].x = -0.744095396034322;
      outPts[15].x = -0.333333333333333;
      outPts[16].x =  0.108864652571010;
      outPts[17].x =  0.468848793470714;
      outPts[18].x = -1.000000000000000;
      outPts[19].x = -0.744095396034321;
      outPts[20].x = -0.364769256536689;
      outPts[21].x =  0.000000000000000;
      outPts[22].x = -1.000000000000000;
      outPts[23].x = -0.769516496193441;
      outPts[24].x = -0.468848793470714;
      outPts[25].x = -1.000000000000000;
      outPts[26].x = -0.830223896278567;
      outPts[27].x = -1.000000000000000;

      outPts[0].y =  -1.000000000000000;
      outPts[1].y =  -1.000000000000000;
      outPts[2].y =  -1.000000000000000;
      outPts[3].y =  -1.000000000000000;
      outPts[4].y =  -1.000000000000000;
      outPts[5].y =  -1.000000000000000;
      outPts[6].y =  -1.000000000000000;
      outPts[7].y =  -0.830223896278567;
      outPts[8].y =  -0.769516496193442;
      outPts[9].y =  -0.744095396034321;
      outPts[10].y = -0.744095396034321;
      outPts[11].y = -0.769516496193441;
      outPts[12].y = -0.830223896278567;
      outPts[13].y = -0.468848793470714;
      outPts[14].y = -0.364769256536689;
      outPts[15].y = -0.333333333333333;
      outPts[16].y = -0.364769256536689;
      outPts[17].y = -0.468848793470714;
      outPts[18].y =  0.000000000000000;
      outPts[19].y =  0.108864652571010;
      outPts[20].y =  0.108864652571010;
      outPts[21].y =  0.000000000000000;
      outPts[22].y =  0.468848793470714;
      outPts[23].y =  0.539032992386883;
      outPts[24].y =  0.468848793470714;
      outPts[25].y =  0.830223896278567;
      outPts[26].y =  0.830223896278567;
      outPts[27].y =  1.000000000000000;
    }
    else if (order == 7) {
      outPts[0].x =  -1.000000000000000;
      outPts[1].x =  -0.871740148509607;
      outPts[2].x =  -0.591700181433142;
      outPts[3].x =  -0.209299217902479;
      outPts[4].x =   0.209299217902479;
      outPts[5].x =   0.591700181433142;
      outPts[6].x =   0.871740148509607;
      outPts[7].x =   1.000000000000000;
      outPts[8].x =  -1.000000000000000;
      outPts[9].x =  -0.823996159427635;
      outPts[10].x = -0.503284033260976;
      outPts[11].x = -0.102529737209732;
      outPts[12].x =  0.304713444152606;
      outPts[13].x =  0.647992318855270;
      outPts[14].x =  0.871740148509606;
      outPts[15].x = -1.000000000000000;
      outPts[16].x = -0.801429410891630;
      outPts[17].x = -0.465073294429312;
      outPts[18].x = -0.069853411141377;
      outPts[19].x =  0.304713444152606;
      outPts[20].x =  0.591700181433142;
      outPts[21].x = -1.000000000000000;
      outPts[22].x = -0.794940525580537;
      outPts[23].x = -0.465073294429312;
      outPts[24].x = -0.102529737209732;
      outPts[25].x =  0.209299217902479;
      outPts[26].x = -1.000000000000000;
      outPts[27].x = -0.801429410891630;
      outPts[28].x = -0.503284033260976;
      outPts[29].x = -0.209299217902479;
      outPts[30].x = -1.000000000000000;
      outPts[31].x = -0.823996159427635;
      outPts[32].x = -0.591700181433142;
      outPts[33].x = -1.000000000000000;
      outPts[34].x = -0.871740148509606;
      outPts[35].x = -1.000000000000000;

      outPts[0].y =  -1.000000000000000;
      outPts[1].y =  -1.000000000000000;
      outPts[2].y =  -1.000000000000000;
      outPts[3].y =  -1.000000000000000;
      outPts[4].y =  -1.000000000000000;
      outPts[5].y =  -1.000000000000000;
      outPts[6].y =  -1.000000000000000;
      outPts[7].y =  -1.000000000000000;
      outPts[8].y =  -0.871740148509606;
      outPts[9].y =  -0.823996159427635;
      outPts[10].y = -0.801429410891630;
      outPts[11].y = -0.794940525580537;
      outPts[12].y = -0.801429410891630;
      outPts[13].y = -0.823996159427635;
      outPts[14].y = -0.871740148509607;
      outPts[15].y = -0.591700181433142;
      outPts[16].y = -0.503284033260976;
      outPts[17].y = -0.465073294429312;
      outPts[18].y = -0.465073294429312;
      outPts[19].y = -0.503284033260976;
      outPts[20].y = -0.591700181433142;
      outPts[21].y = -0.209299217902479;
      outPts[22].y = -0.102529737209732;
      outPts[23].y = -0.069853411141377;
      outPts[24].y = -0.102529737209731;
      outPts[25].y = -0.209299217902479;
      outPts[26].y =  0.209299217902479;
      outPts[27].y =  0.304713444152606;
      outPts[28].y =  0.304713444152606;
      outPts[29].y =  0.209299217902479;
      outPts[30].y =  0.591700181433142;
      outPts[31].y =  0.647992318855270;
      outPts[32].y =  0.591700181433142;
      outPts[33].y =  0.871740148509606;
      outPts[34].y =  0.871740148509606;
      outPts[35].y =  1.000000000000000;
    }
    else if (order == 8) {
      outPts[0].x =  -1.000000000000000;
      outPts[1].x =  -0.899757995411460;
      outPts[2].x =  -0.677186279510738;
      outPts[3].x =  -0.363117463826178;
      outPts[4].x =   0.000000000000000;
      outPts[5].x =   0.363117463826178;
      outPts[6].x =   0.677186279510738;
      outPts[7].x =   0.899757995411460;
      outPts[8].x =   1.000000000000000;
      outPts[9].x =  -1.000000000000000;
      outPts[10].x = -0.861819860276539;
      outPts[11].x = -0.602543502716998;
      outPts[12].x = -0.264723799928576;
      outPts[13].x =  0.101891561713898;
      outPts[14].x =  0.446427110295541;
      outPts[15].x =  0.723639720553079;
      outPts[16].x =  0.899757995411460;
      outPts[17].x = -1.000000000000000;
      outPts[18].x = -0.843883607578542;
      outPts[19].x = -0.566537963017939;
      outPts[20].x = -0.221988381746101;
      outPts[21].x =  0.133075926035877;
      outPts[22].x =  0.446427110295541;
      outPts[23].x =  0.677186279510738;
      outPts[24].x = -1.000000000000000;
      outPts[25].x = -0.837167761785322;
      outPts[26].x = -0.556023236507798;
      outPts[27].x = -0.221988381746101;
      outPts[28].x =  0.101891561713898;
      outPts[29].x =  0.363117463826178;
      outPts[30].x = -1.000000000000000;
      outPts[31].x = -0.837167761785322;
      outPts[32].x = -0.566537963017939;
      outPts[33].x = -0.264723799928576;
      outPts[34].x =  0.000000000000000;
      outPts[35].x = -1.000000000000000;
      outPts[36].x = -0.843883607578542;
      outPts[37].x = -0.602543502716998;
      outPts[38].x = -0.363117463826178;
      outPts[39].x = -1.000000000000000;
      outPts[40].x = -0.861819860276539;
      outPts[41].x = -0.677186279510738;
      outPts[42].x = -1.000000000000000;
      outPts[43].x = -0.899757995411460;
      outPts[44].x = -1.000000000000000;

      outPts[0].y =  -1.000000000000000;
      outPts[1].y =  -1.000000000000000;
      outPts[2].y =  -1.000000000000000;
      outPts[3].y =  -1.000000000000000;
      outPts[4].y =  -1.000000000000000;
      outPts[5].y =  -1.000000000000000;
      outPts[6].y =  -1.000000000000000;
      outPts[7].y =  -1.000000000000000;
      outPts[8].y =  -1.000000000000000;
      outPts[9].y =  -0.899757995411460;
      outPts[10].y = -0.861819860276539;
      outPts[11].y = -0.843883607578542;
      outPts[12].y = -0.837167761785322;
      outPts[13].y = -0.837167761785322;
      outPts[14].y = -0.843883607578542;
      outPts[15].y = -0.861819860276539;
      outPts[16].y = -0.899757995411460;
      outPts[17].y = -0.677186279510738;
      outPts[18].y = -0.602543502716999;
      outPts[19].y = -0.566537963017939;
      outPts[20].y = -0.556023236507798;
      outPts[21].y = -0.566537963017938;
      outPts[22].y = -0.602543502716999;
      outPts[23].y = -0.677186279510738;
      outPts[24].y = -0.363117463826178;
      outPts[25].y = -0.264723799928576;
      outPts[26].y = -0.221988381746101;
      outPts[27].y = -0.221988381746101;
      outPts[28].y = -0.264723799928576;
      outPts[29].y = -0.363117463826178;
      outPts[30].y =  0.000000000000000;
      outPts[31].y =  0.101891561713898;
      outPts[32].y =  0.133075926035877;
      outPts[33].y =  0.101891561713898;
      outPts[34].y =  0.000000000000000;
      outPts[35].y =  0.363117463826178;
      outPts[36].y =  0.446427110295541;
      outPts[37].y =  0.446427110295541;
      outPts[38].y =  0.363117463826178;
      outPts[39].y =  0.677186279510738;
      outPts[40].y =  0.723639720553079;
      outPts[41].y =  0.677186279510738;
      outPts[42].y =  0.899757995411460;
      outPts[43].y =  0.899757995411460;
      outPts[44].y =  1.000000000000000;
    }
  }
  else if (eType == QUAD) {
    // Tensor-product element
    vector<double> spts1D = getPts1D(params->sptsTypeQuad,order);
    outPts.resize((order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        outPts[j+i*(order+1)].x = spts1D[j];
        outPts[j+i*(order+1)].y = spts1D[i];
      }
    }
  }
  else if (eType == HEX) {
    // Tensor-product element
    vector<double> spts1D = getPts1D(params->sptsTypeQuad,order);
    outPts.resize((order+1)*(order+1)*(order+1));
    for (int k=0; k<order+1; k++) {
      for (int j=0; j<order+1; j++) {
        for (int i=0; i<order+1; i++) {
          outPts[i+(order+1)*(j+(order+1)*k)].x = spts1D[i];
          outPts[i+(order+1)*(j+(order+1)*k)].y = spts1D[j];
          outPts[i+(order+1)*(j+(order+1)*k)].z = spts1D[k];
        }
      }
    }
  }

  return outPts;
}

vector<point> geo::getLocFpts(int eType, int order)
{
  vector<point> outPts;
  vector<double> pts1D;

  if (eType == TRI) {
    outPts.resize(3*(order+1));
    pts1D = getPts1D(params->sptsTypeTri,order);
    for (int i=0; i<order+1; i++) {
      // Face 0
      outPts[i].x = pts1D[i];
      outPts[i].y = -1.;
      // Face 1
      outPts[i+order+1].x = pts1D[order-i];
      outPts[i+order+1].y = pts1D[i];
      // Face 2
      outPts[i+2*(order+1)].x = -1.;
      outPts[i+2*(order+1)].y = pts1D[order-i];
    }
  }
  else if (eType == QUAD) {
    outPts.resize(4*(order+1));
    pts1D = getPts1D(params->sptsTypeQuad,order);
    for (int i=0; i<order+1; i++) {
      // Face 0
      outPts[i].x = pts1D[i];
      outPts[i].y = -1.;
      // Face 1
      outPts[i+order+1].x = 1.;
      outPts[i+order+1].y = pts1D[i];
      // Face 2
      outPts[i+2*(order+1)].x = pts1D[order-i];
      outPts[i+2*(order+1)].y = 1.;
      // Face 3
      outPts[i+3*(order+1)].x = -1.;
      outPts[i+3*(order+1)].y = pts1D[order-i];
    }
  }
  else if (eType == HEX) {
    int P12 = (order+1)*(order+1);
    outPts.resize(6*P12);
    pts1D = getPts1D(params->sptsTypeQuad,order);
    // Flux points are ordered such that, as seen from inside the
    // element, the id's increase btm-left->top-right fashion on
    // each face, starting with lowest dimension first ('x' or 'y')
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        int ind = i+j*(order+1);
        // Face 0 - bottom
        outPts[ind].x = pts1D[i];
        outPts[ind].y = pts1D[j];
        outPts[ind].z = -1.;
        // Face 1 - top
        outPts[ind+P12].x = pts1D[order-i];
        outPts[ind+P12].y = pts1D[j];
        outPts[ind+P12].z = 1.;
        // Face 2 - left
        outPts[ind+2*P12].x = -1;
        outPts[ind+2*P12].y = pts1D[i];
        outPts[ind+2*P12].z = pts1D[j];
        // Face 3 - right
        outPts[ind+3*P12].x = 1;
        outPts[ind+3*P12].y = pts1D[order-i];
        outPts[ind+3*P12].z = pts1D[j];
        // Face 4 - front
        outPts[ind+4*P12].x = pts1D[order-i];
        outPts[ind+4*P12].y = -1;
        outPts[ind+4*P12].z = pts1D[j];
        // Face 5 - back
        outPts[ind+5*P12].x = pts1D[i];
        outPts[ind+5*P12].y = 1;
        outPts[ind+5*P12].z = pts1D[j];
      }
    }
  }

  return outPts;
}

vector<double> geo::getPts1D(string ptsType, int order)
{
  vector<double> outPts(order+1);

  if (!ptsType.compare("Legendre")) { // Gauss-Legendre
    if (order == 0) {
      outPts[0] =  0.0;
    }
    else if(order == 1) {
      outPts[0] = -0.577350269189626;
      outPts[1] =  0.577350269189626;
    }
    else if(order == 2) {
      outPts[0] = -0.774596669241483;
      outPts[1] =  0.000000000000000;
      outPts[2] =  0.774596669241483;
    }
    else if(order == 3) {
      outPts[0] = -0.861136311594053;
      outPts[1] = -0.339981043584856;
      outPts[2] =  0.339981043584856;
      outPts[3] =  0.861136311594053;
    }
    else if(order == 4) {
      outPts[0] = -0.906179845938664;
      outPts[1] = -0.538469310105683;
      outPts[2] =  0.000000000000000;
      outPts[3] =  0.538469310105683;
      outPts[4] =  0.906179845938664;
    }
    else if(order == 5) {
      outPts[0] = -0.932469514203152;
      outPts[1] = -0.661209386466264;
      outPts[2] = -0.238619186083197;
      outPts[3] =  0.238619186083197;
      outPts[4] =  0.661209386466265;
      outPts[5] =  0.932469514203152;
    }
    else if(order == 6) {
      outPts[0] = -0.949107912342758;
      outPts[1] = -0.741531185599395;
      outPts[2] = -0.405845151377397;
      outPts[3] =  0.000000000000000;
      outPts[4] =  0.405845151377397;
      outPts[5] =  0.741531185599394;
      outPts[6] =  0.949107912342758;
    }
    else if(order == 7) {
      outPts[0] = -0.960289856497536;
      outPts[1] = -0.796666477413627;
      outPts[2] = -0.525532409916329;
      outPts[3] = -0.183434642495650;
      outPts[4] =  0.183434642495650;
      outPts[5] =  0.525532409916329;
      outPts[6] =  0.796666477413627;
      outPts[7] =  0.960289856497536;
    }
    else if(order == 8) {
      outPts[0] = -0.968160239507626;
      outPts[1] = -0.836031107326636;
      outPts[2] = -0.613371432700591;
      outPts[3] = -0.324253423403809;
      outPts[4] =  0.000000000000000;
      outPts[5] =  0.324253423403809;
      outPts[6] =  0.613371432700591;
      outPts[7] =  0.836031107326636;
      outPts[8] =  0.968160239507626;
    }
    else if(order == 9) {
      outPts[0] = -0.973906528517172;
      outPts[1] = -0.865063366688985;
      outPts[2] = -0.679409568299024;
      outPts[3] = -0.433395394129247;
      outPts[4] = -0.148874338981631;
      outPts[5] =  0.148874338981632;
      outPts[6] =  0.433395394129247;
      outPts[7] =  0.679409568299024;
      outPts[8] =  0.865063366688984;
      outPts[9] =  0.973906528517171;
    }
    else if(order == 10) {
      outPts[0] = -0.978228658146057;
      outPts[1] = -0.887062599768095;
      outPts[2] = -0.730152005574049;
      outPts[3] = -0.519096129206812;
      outPts[4] = -0.269543155952345;
      outPts[5] =  0.000000000000000;
      outPts[6] =  0.269543155952345;
      outPts[7] =  0.519096129206812;
      outPts[8] =  0.730152005574049;
      outPts[9] =  0.887062599768095;
      outPts[10]=  0.978228658146057;
    }
    else {
      stringstream ss; ss << order;
      string errMsg = "Gauss-Ledgendre point locations for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }
  else if (!ptsType.compare("Lobatto")) { // Gauss-Lobatto
    if (order == 0) {
        outPts[0] =  0.0;
    }
    else if(order == 1) {
      outPts[0] = -1.000000000000000;
      outPts[1] =  1.000000000000000;
    }
    else if(order == 2) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.000000000000000;
      outPts[2] =  1.000000000000000;
    }
    else if(order == 3) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.447213595499958;
      outPts[2] =  0.447213595499958;
      outPts[3] =  1.000000000000000;
    }
    else if(order == 4) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.654653670707977;
      outPts[2] = -0.000000000000000;
      outPts[3] =  0.654653670707977;
      outPts[4] =  1.000000000000000;
    }
    else if(order == 5) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.765055323929465;
      outPts[2] = -0.285231516480645;
      outPts[3] =  0.285231516480645;
      outPts[4] =  0.765055323929465;
      outPts[5] =  1.000000000000000;
    }
    else if(order == 6) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.830223896278567;
      outPts[2] = -0.468848793470714;
      outPts[3] =  0.000000000000000;
      outPts[4] =  0.468848793470714;
      outPts[5] =  0.830223896278567;
      outPts[6] =  1.000000000000000;
    }
    else if(order == 7) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.871740148509607;
      outPts[2] = -0.591700181433142;
      outPts[3] = -0.209299217902479;
      outPts[4] =  0.209299217902479;
      outPts[5] =  0.591700181433142;
      outPts[6] =  0.871740148509607;
      outPts[7] =  1.000000000000000;
    }
    else if(order == 8) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.899757995411460;
      outPts[2] = -0.677186279510738;
      outPts[3] = -0.363117463826178;
      outPts[4] = -0.000000000000000;
      outPts[5] =  0.363117463826178;
      outPts[6] =  0.677186279510738;
      outPts[7] =  0.899757995411460;
      outPts[8] =  1.000000000000000;
    }
    else if(order == 9) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.919533908166459;
      outPts[2] = -0.738773865105505;
      outPts[3] = -0.477924949810445;
      outPts[4] = -0.165278957666387;
      outPts[5] =  0.165278957666387;
      outPts[6] =  0.477924949810444;
      outPts[7] =  0.738773865105505;
      outPts[8] =  0.919533908166459;
      outPts[9] =  1.000000000000000;
    }
    else if(order == 10) {
      outPts[0] = -1.000000000000000;
      outPts[1] = -0.934001430408059;
      outPts[2] = -0.784483473663144;
      outPts[3] = -0.565235326996205;
      outPts[4] = -0.295758135586939;
      outPts[5] =  0.000000000000000;
      outPts[6] =  0.295758135586939;
      outPts[7] =  0.565235326996205;
      outPts[8] =  0.784483473663144;
      outPts[9] =  0.934001430408059;
      outPts[10] =  1.000000000000000;
    }
    else {
      stringstream ss;  ss << order;
      string errMsg = "Gauss-Lobatto point locations for order " + ss.str() + " not implemented.";
      FatalError(errMsg.c_str());
    }
  }

  return outPts;
}

vector<double> geo::getQptWeights(int order)
{
  // Tensor-product elements
  vector<double> qwts1D = getQptWeights1D(order);
  vector<double> outWts;
  if (nDims == 2) {
    outWts.resize((order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        outWts[j+i*(order+1)] = qwts1D[i]*qwts1D[j];
      }
    }
  }
  else if (nDims == 3) {
    outWts.resize((order+1)*(order+1)*(order+1));
    for (int i=0; i<order+1; i++) {
      for (int j=0; j<order+1; j++) {
        for (int k=0; k<order+1; k++) {
          outWts[k+(order+1)*(j+i*(order+1))] = qwts1D[i]*qwts1D[j]*qwts1D[k];
        }
      }
    }
  }

  return outWts;
}


vector<double> geo::getQptWeights1D(int order)
{
  // Order here refers to the order of a polynomial fit through
  // the Gauss points, not the order of accuracy of integration
  // using the same number of points

  vector<double> outWts(order+1);

  if (order == 0) {
    outWts[0] =  2.0;
  }
  else if(order == 1) {
    outWts[0] = 1.0;
    outWts[1] = 1.0;
  }
  else if(order == 2) {
    outWts[0] = 0.5555555555555556;
    outWts[1] = 0.8888888888888888;
    outWts[2] = 0.5555555555555556;
  }
  else if(order == 3) {
    outWts[0] = 0.3478548451374538;
    outWts[1] = 0.6521451548625461;
    outWts[2] = 0.6521451548625461;
    outWts[3] = 0.3478548451374538;
  }
  else if(order == 4) {
    outWts[0] = 0.2369268850561891;
    outWts[1] = 0.4786286704993665;
    outWts[2] = 0.000000000000000;
    outWts[3] = 0.4786286704993665;
    outWts[4] = 0.2369268850561891;
  }
  else {
    stringstream ss; ss << order;
    string errMsg = "Gauss quadrature weights for order " + ss.str() + " not implemented.";
    FatalError(errMsg.c_str());
  }

  return outWts;
}
