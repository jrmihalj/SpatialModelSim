int Dispersal_Func(float *N, float *Z, float *B, float *oldN, float *oldZ, float *oldB,
                   int Lat_Size, int Cells, double disperse, gsl_rng *r2){
  
  
  int c;
  int Num_neighbors[Cells];
  int neighbor;
  int Row, Col; //For looping
  
  float *Subtract_N, *Subtract_Z;
  float *Add_N, *Add_Z;
  Subtract_N = vector(1,Cells);
  Subtract_Z = vector(1,Cells);
  Add_N = vector(1,Cells);
  Add_Z = vector(1,Cells);
  
  
  Row=1; Col=1;
  //Figure out subtractions:
  //How many leaving the cell?
  for(c=1;c<=Cells;c++){
    
    // = Current(old)Density * %Disperse (up to disperse)
    Subtract_N[c] = oldN[c] * disperse * gsl_rng_uniform(r2) / 2;
    Subtract_Z[c] = oldZ[c] * disperse * gsl_rng_uniform(r2) / 2;
    
    if( Row==1 || Row==Lat_Size ){ //If the top or bottom row
      Num_neighbors[c] = 5;
    }else{
      Num_neighbors[c] = 8;
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
    Add_N[c] = 0; Add_Z[c] = 0; //Set all equal to zero to start
  }
  
  for(c=1;c<=Cells;c++){
    
    //Right neighbors//:
    //Immediate right and bottom right
    if( Col==Lat_Size ){
      
      neighbor = c - Lat_Size - 1; //adjacent right
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      
      if(Row != Lat_Size){
        neighbor = c + 1;
        Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      }
 
    }else{
      neighbor = c + 1;
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
      
      neighbor = c + Lat_Size + 1;
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
    }
    
    //Left neighbors:
    if( Col==1 ){
      neighbor = c + (Lat_Size-1);
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
    }else{
      neighbor = c - 1;
      Add_N[c] += Subtract_N[neighbor] / Num_neighbors[neighbor];
    }
    
    
    
    
      
    //check for end of row:
    if( c % Lat_Size == 0){
      Row++;
      Col=0;
    }
    
    //Progress one column:
    Col++;
      
  }
  
  
  free_vector(Subtract_N,1,Cells);
  
  
}