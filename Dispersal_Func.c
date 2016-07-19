int Dispersal_Func(float *oldN, float *oldZ, float *oldB,
                   int Lat_Size, int Cells, double disperse, gsl_rng *r2, int y){
  
  FILE *fp4;
  
  int c;
  int Num_neighbors[Cells];
  int neighbor;
  int Row, Col; //For looping
  
  float *Subtract_N, *Subtract_Z;
  float *Add_N, *Add_Z;
  Subtract_N = vector(1,Cells);
  //Subtract_Z = vector(1,Cells);
  Add_N = vector(1,Cells);
  Add_Z = vector(1,Cells);
  
  
  Row=1; Col=1;
  //Figure out subtractions:
  //How many leaving the cell?
  for(c=1;c<=Cells;c++){
    
    // = Current(old) Density * %Disperse (up to disperse) / 2 (only males disperse)
    Subtract_N[c] = oldN[c] * disperse * gsl_rng_uniform(r2) / 2;
    //Subtract_Z[c] = oldZ[c] * disperse * gsl_rng_uniform(r2) / 2;
    
    
    
    if( Row==1 || Row==Lat_Size ){ 
      
      
      if( Col==1 || Col==Lat_Size ){ //If a corner cell
        Num_neighbors[c] = 3;
      }else{ //If top or bottom edge
        Num_neighbors[c] = 5;
      }
      
    }else{
      
      if( Col==1 || Col==Lat_Size ){ //If side edge
        Num_neighbors[c] = 5;
      }else{ //If interior
        Num_neighbors[c] = 8;
      }
      
    }
    
    //check for end of row:
    if( c % Lat_Size == 0){
      Row++;
      Col=0;
    }
    
    //Progress one column:
    Col++;
    
  }
  
  /**************************************************/
  /**************************************************/
  
  //Figure out additions:
  
  for(c=1;c<=Cells;c++){
    Add_N[c] = 0; //Add_Z[c] = 0; //Set all equal to zero to start
  }
  
  Row=1; Col=1;
  
  for(c=1;c<=Cells;c++){
    
    //Right neighbors//:
    if( Col != Lat_Size){
      
      //adjacent right
      neighbor = c + 1; 
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      
      //upper right
      if( Row != 1 ){
        neighbor = c - Lat_Size + 1; 
        Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
        //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      }
      
      //lower right
      if( Row != Lat_Size ){
        neighbor = c + Lat_Size + 1; 
        Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
        //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      }
      
    }
    
    //Left neighbors//:
    if( Col != 1 ){
      
      //adjacent left
      neighbor = c - 1; 
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      
      //upper left
      if( Row != 1 ){
        neighbor = c - Lat_Size - 1; 
        Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
        //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      }
      
      //lower left
      if( Row != Lat_Size ){
        neighbor = c + Lat_Size - 1; 
        Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
        //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
      }
      
    }
    
    //Top neighbors//:
    if(Row != 1){
      neighbor = c - Lat_Size; 
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
    }
    
    //Bottom neighbors//:
    if(Row != Lat_Size){
      neighbor = c + Lat_Size; 
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      //Add_Z[c] += Subtract_Z[neighbor] / Num_neighbors[neighbor];
    }
    
    
    // UPDATE ROW AND COLUMN ID //
    
    
    //check for end of row:
    if( c % Lat_Size == 0){
      Row++;
      Col=0;
    }
    
    //Progress one column:
    Col++;
      
  }
  
  /**************************************************/
  /**************************************************/
  
  //Add these up:
  //Also print:
  
  if(Debug==1){
    printf("%s\n",
           "Starting dispersal update...");
  }
  
  for(c=1;c<=Cells;c++){
    
    if(Debug==1){
      printf("%s\t %d\t %d\t %f\n",
             "before", y, c, oldN[c]);
    }
    
    oldN[c] = oldN[c] + Add_N[c] - Subtract_N[c];
    //oldZ[c] = oldZ[c] + Add_Z[c] - Subtract_Z[c];
    
    if(Debug==1){
      printf("%s\t %d\t %d\t %f\n",
             "after", y, c, oldN[c]);
    }
    
  }
  
  //Still need to do something with Biomass...
  
  free_vector(Subtract_N,1,Cells);
  free_vector(Subtract_Z,1,Cells);
  free_vector(Add_N,1,Cells);
  free_vector(Add_Z,1,Cells);
  
  
}