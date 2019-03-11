#include "cn.h"
int ymw16par(struct Warp_Sun *t0, struct Thick *t1, struct Thin *t2, struct Spiral *t3, struct GC *t4, struct Gum *t5,struct LB *t6,  struct LI *t7, struct FB *t8, struct  LMC *t9, struct Dora *t10, struct SMC *t11, char *dirname){

  //FILE *fptr =NULL;
  //char key[40], *cstr, filen[64];
  //double value;
  //double a[5];
  //int i;
  //size_t size=100;

  //cstr = (char *)malloc(sizeof(char)*size);
  
  //strcpy(filen,dirname);
  //strcat(filen,"ymw16par.txt");
  //fptr = fopen(filen,"r");

  //if(fptr == NULL)
  //{
  //  printf("File %s open error\n",filen);
  //  return 0;
  //}

  t0->Gamma_w = 0.14;
  t0->z_Sun = 6.0;

  //printf("TEST for constant: %8.3f %8.3f", t0->Gamma_w, t0->z_Sun);


  t1->Ad = 2500.0;
  t1->Bd = 15000.0;
  t1->n1 = 0.01132;
  t1->H1 = 1673.0;

  t2->A2 = 1200.0;
  t2->B2 = 4000.0;
  t2->n2 = 0.404;
  t2->K2 = 1.54;

  t3->B2s = 4000.0;
  t3->Ka  = 5.015;
  t3->narm[0] = 0.135;
  t3->narm[1] = 0.129;
  t3->narm[2] = 0.103;
  t3->narm[3] = 0.116;
  t3->narm[4] = 0.0057;
  t3->warm[0] = 300.0;
  t3->warm[1] = 500.0;
  t3->warm[2] = 300.0;
  t3->warm[3] = 500.0;
  t3->warm[4] = 300.0;
  t3->Aa = 11680.0;
  t3->ncn = 2.4;
  t3->wcn = 8.2;
  t3->thetacn = 109.0;
  t3->nsg = 0.626;
  t3->wsg = 20.0;
  t3->thetasg = 75.8;

  t4->ngc = 6.2;
  t4->Agc = 160.0;
  t4->Hgc = 35.0;
  
  t5->ngn = 1.84;
  t5->Wgn = 15.1;
  t5->Agn = 125.8;
  t5->Kgn = 1.4;
  
  t6->J_LB = 0.480;
  t6->nlb1 = 1.094;
  t6->detlb1 = 28.4;
  t6->wlb1 = 14.2;
  t6->hlb1 = 112.9;
  t6->thetalb1 = 195.4;
  t6->nlb2 = 2.33;
  t6->detlb2 = 14.7;
  t6->wlb2 = 15.6;
  t6->hlb2 = 43.6;
  t6->thetalb2 = 278.2;

  t7->nLI = 1.907;
  t7->RLI = 80.0;
  t7->WLI = 15.0;
  t7->detthetaLI = 30.0;
  t7->thetaLI = 40.0;

  t8->J_FB = 1.0;

  t9->nlmc = 0.066;

  t10->n30D = 0.32;

  t11->nsmc = 0.045;

 /*
  while(!feof(fptr))
    {
      if( fscanf(fptr,"%s",key) == 1){
	if (key[0]=='#'){
	  getline(&cstr,&size,fptr);
	  //        printf("%s \n",cstr);
	}
	
	//Warp and Sun
    else if(strcmp(key,"Gamma_w") == 0){
	  fscanf(fptr,"%lf",&((*t0).Gamma_w ) );
	}    
	else if(strcmp(key,"z_Sun") == 0){
	  fscanf(fptr,"%lf",&((*t0).z_Sun ) );
	}

	//thick disk
	else if(strcmp(key,"Ad") == 0){
	  fscanf(fptr,"%lf",&((*t1).Ad ) );
	}    
	else if(strcmp(key,"Bd") == 0){
	  fscanf(fptr,"%lf",&((*t1).Bd ) );
	}
	
	else if(strcmp(key,"n1") == 0){
	  fscanf(fptr,"%lf",&((*t1).n1 ) );
	}
	else if(strcmp(key,"H1") == 0){
	  fscanf(fptr,"%lf",&((*t1).H1 ) );
	} 
    //thin disk
	else if(strcmp(key,"A2") == 0){
	  fscanf(fptr,"%lf",&((*t2).A2 ) );
	}
	else if(strcmp(key,"B2") == 0){
	  fscanf(fptr,"%lf",&((*t2).B2 ) );
	}
	
	else if(strcmp(key,"n2") == 0){
	  fscanf(fptr,"%lf",&((*t2).n2 ) );
	}
	else if(strcmp(key,"K2") == 0){
	  fscanf(fptr,"%lf",&((*t2).K2 ) );
	}
	
	//spiral arm
	else if(strcmp(key,"B2s") == 0){
	  fscanf(fptr,"%lf",&((*t3).B2s ) );    
	}
	else if(strcmp(key,"Aa") == 0){
	  fscanf(fptr,"%lf",&((*t3).Aa ) );    
	}
	else if(strcmp(key,"ncn") == 0){
	  fscanf(fptr,"%lf",&((*t3).ncn) );        
	}  
	else if(strcmp(key,"wcn") == 0){
	  fscanf(fptr,"%lf",&((*t3).wcn ) );
	}
	else if(strcmp(key,"thetacn") == 0){
	  fscanf(fptr,"%lf",&((*t3).thetacn ) );
	}
	else if(strcmp(key,"nsg") == 0){
	  fscanf(fptr,"%lf",&((*t3).nsg ) );
	}
	else if(strcmp(key,"wsg") == 0){
	  fscanf(fptr,"%lf",&((*t3).wsg ) );
	}
	else if(strcmp(key,"thetasg") == 0){
	  fscanf(fptr,"%lf",&((*t3).thetasg ) );
	}
	else if(strcmp(key,"Ka") == 0){
	  fscanf(fptr,"%lf",&((*t3).Ka ) );
	}
	else if(strcmp(key,"Ele_arm") == 0){
	  fscanf(fptr,"%lf %lf %lf %lf %lf",&((*t3).narm[0]),&((*t3).narm[1]),&((*t3).narm[2]),&((*t3).narm[3]),&((*t3).narm[4])  );
	}
	else if(strcmp(key,"Wid_arm") == 0){
        fscanf(fptr,"%lf %lf %lf %lf %lf",&((*t3).warm[0]),&((*t3).warm[1]),&((*t3).warm[2]),&((*t3).warm[3]),&((*t3).warm[4])  );
	}
      
     
  
	//Galactic center
	else if(strcmp(key,"Agc") == 0){
	  fscanf(fptr,"%lf",&((*t4).Agc ) ); 
	}
	else if(strcmp(key,"Hgc") == 0){
	  fscanf(fptr,"%lf",&((*t4).Hgc ) );                   
	}
	else if(strcmp(key,"ngc") == 0){
	  fscanf(fptr,"%lf",&((*t4).ngc ) ); 
	}

	//Gum nebula
	else if(strcmp(key,"Kgn") == 0){
	  fscanf(fptr,"%lf",&((*t5).Kgn ) );
	}
	else if(strcmp(key,"ngn") == 0){
	  fscanf(fptr,"%lf",&((*t5).ngn ) );
	}
	else if(strcmp(key,"Wgn") == 0){
	  fscanf(fptr,"%lf",&((*t5).Wgn ) );
	}
	else if(strcmp(key,"Agn") == 0){
	  fscanf(fptr,"%lf",&((*t5).Agn ) );
	}
	//Local Bubble 
	else if(strcmp(key,"J_LB") == 0){
	  fscanf(fptr,"%lf",&((*t6).J_LB ) );
	}
	else if(strcmp(key,"nlb1") == 0){
	  fscanf(fptr,"%lf",&((*t6).nlb1 ) );
	}
	else if(strcmp(key,"detlb1") == 0){
	  fscanf(fptr,"%lf",&((*t6).detlb1 ) );
	}
	else if(strcmp(key,"wlb1") == 0){
	  fscanf(fptr,"%lf",&((*t6).wlb1 ) );
	}
	else if(strcmp(key,"hlb1") == 0){
	  fscanf(fptr,"%lf",&((*t6).hlb1 ) );
	}
	else if(strcmp(key,"thetalb1") == 0){
	  fscanf(fptr,"%lf",&((*t6).thetalb1 ) );
	}
	else if(strcmp(key,"nlb2") == 0){
	  fscanf(fptr,"%lf",&((*t6).nlb2 ) );
	}
	else if(strcmp(key,"detlb2") == 0){
	  fscanf(fptr,"%lf",&((*t6).detlb2 ) );
	}
	else if(strcmp(key,"wlb2") == 0){
	  fscanf(fptr,"%lf",&((*t6).wlb2 ) );
	}
	else if(strcmp(key,"hlb2") == 0){
	  fscanf(fptr,"%lf",&((*t6).hlb2 ) );
	}
	else if(strcmp(key,"thetalb2") == 0){
	  fscanf(fptr,"%lf",&((*t6).thetalb2 ) );
	}
	
	//Loop I
	else if(strcmp(key,"nLI") == 0){
	  fscanf(fptr,"%lf",&((*t7).nLI ) );
	}
	else if(strcmp(key,"RLI") == 0){
	  fscanf(fptr,"%lf",&((*t7).RLI ) );
	}
	else if(strcmp(key,"WLI") == 0){
	  fscanf(fptr,"%lf",&((*t7).WLI ) );
	}
	else if(strcmp(key,"detthetaLI") == 0){
	  fscanf(fptr,"%lf",&((*t7).detthetaLI ) );
	}
	else if(strcmp(key,"thetaLI") == 0){
	  fscanf(fptr,"%lf",&((*t7).thetaLI ) );
	}

	//Fermi Bubble
	else if(strcmp(key,"J_FB") == 0){
	  fscanf(fptr,"%lf",&((*t8).J_FB ) );
	}
	
	//LMC
	else if(strcmp(key,"nlmc") == 0){
	  fscanf(fptr,"%lf",&((*t9).nlmc ) );
	}
	
	//30 Doradus
	else if(strcmp(key,"n30D") == 0){
	  fscanf(fptr,"%lf",&((*t10).n30D ) );
	}
	//SMC
	else if(strcmp(key,"nsmc") == 0){
	  fscanf(fptr,"%lf",&((*t11).nsmc ) );
	}
      }
    }
 
printf("test for reading all the parameters: %8.3f, %8.3f, %8.3f", t0->Gamma_w, t1->Ad, t10->n30D);
*/
 
//  free(cstr);
}

