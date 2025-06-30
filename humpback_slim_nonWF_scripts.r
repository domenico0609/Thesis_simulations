library(slimr)

##Slim script for the simulation of the nuclear and mitochondrial DNA. Nm=0
slim_script(
  slim_block(initialize(), 
             {
               ##initializing a nucleotide based non-WF simulation
               initializeSLiMModelType("nonWF");	
               initializeSLiMOptions(nucleotideBased=T);
               
               ##reading a Fasta file as Ancestral string for the individuals
               defineConstant("L", initializeAncestralNucleotides("C:/Users/Domenico/Documents/humpy/Ancestral_Fasta_test.fasta"));
               
               ##Sexual reproduction call
               initializeSex("A");	
               
               ##vector of carrying capacities for each subpopulation and a vector for the benefit of the specific feeding ground of the subpopulation
               defineConstant("K", c(1000.0, 1000.0, 1000.0, 1000.0));	
               
               ##number of populations
               defineConstant("P", 4);
               
               ##Number of individuals per subpopulation at the start of the simulation 
               defineConstant("N", 10);	
               
               ##Age of maturity and maximum age of the population are defined here (will be changed to a sex-based age of maturity later on. This would require additionial changed to the mortality age table
               defineConstant("Age_mat", 5);
               defineConstant("Age_max", 80);
               
               ##to counteract the early boost in fitness from the carrying capacity formula a low-density value is defined so that at low density there is no additional advantage for an individual. This approaches a more natural increase in individuals in contrast to the early boost from the carrying capacity formula.
               defineConstant("Low_density_benefit", asFloat(1.0));
               
               #Sex based survivability
               defineConstant("Sur_male", 0.98);
               defineConstant("Sur_female", 0.95);
               
               Mort = float(Age_max+1);
               Mort[(Age_max-10):Age_max] = c(0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);
               Mort[0:Age_mat] = c(0.05, 0.05, 0.025, 0.025, 0.025, 0.025);
               
               defineConstant ("Mort_age", Mort);
               
               for (i in 1:10) initializeMutationTypeNuc(i, 0.5, "f", 0.0);
               initializeMutationType("m11", 0.5, "f", 0.0);
               
               m1.color = "purple";
               m2.color = "pink";
               m3.color = "blue";
               m4.color = "red";
               m5.color = "yellow";
               m6.color = "green";
               m7.color = "blue";
               m8.color = "red";
               m9.color = "yellow";
               m10.color = "green";
               m11.color = "grey";
               
               initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(1.12e-8));	##neutral nuclear DNA region
               initializeGenomicElementType("g2", m2, 1.0, mmJukesCantor(4.3e-6));	##neutral mitochondrial DNA region 
               initializeGenomicElement(g1, 0, asInteger(L-22));
               initializeGenomicElement(g2, asInteger(L-21), asInteger(L-1));
               initializeRecombinationRate(c(1e-8, 0), c(asInteger(L-22), asInteger(L-1)));	
               
               g1.color = "blue";
               g2.color = "pink";	
               
             }
  ),
  
  slim_block(reproduction(), 
             
             {	
               if (individual.sex == "F" & individual.age >= Age_mat & individual.tag == 0) {
                 
                 matePop = sample(sim.subpopulations, 1);
                 mate= matePop.sampleIndividuals(1, sex = "M", minAge = Age_mat);
                 
                 ##Checks if the above sample actually contains a male
                 if (size(mate) > 0)
                 { 
                   child = subpop.addCrossed(individual, mate);
                   child.tag = 0;
                 }
                 
                 individual.tag = rdunif(1, min=1, max=2);
               }
               
               ##if the individual has a tag higher than 0 it is lowered by 1 and once the tag reaches zero again the female can reproduce again
               else if (individual.tag > 0) individual.tag = individual.tag - 1;
             }
  ),
  
  slim_block(1, early(),
             {
               for (i in 1:P) sim.addSubpop(i, N);
               
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1N.vcf", m3);
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1M.vcf", m7);
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2N.vcf", m4);	
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2M.vcf", m8);
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3N.vcf", m5);	
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3M.vcf", m9);
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4N.vcf", m6);	
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4M.vcf", m10);
               
               for (subpop in sim.subpopulations) 
               {
                 ind = subpop.individuals;
                 ind.age = rdunif(N, min=0, max=Age_max-10);
                 ind.tag = 0;
               }		
               
               gen1 = sim.subpopulations.individuals;
               gen1.genome1.addNewMutation(m11, 0.0, L-1);
               
             }
  ),
  
  slim_block(early(),
             {
               o = 0;
               
               ##life table mortality p1
               for (subpop in sim.subpopulations)
               {
                 if (subpop.individualCount > 0) {
                   inds = subpop.individuals;
                   ages = inds.age;
                   age_mortality = Mort_age[ages];
                   age_survival = 1 - age_mortality;
                   
                   ##fitness is calculated for all individuals based on their age and feeding ground
                   inds.fitnessScaling = age_survival;
                   
                   ##the fitnessScaling is altered by a survival factor based on the sex of the individual
                   females = 0;
                   males = 0;		
                   
                   for (ind in inds)
                     if (ind.sex =="F"){
                       ind.fitnessScaling = ind.fitnessScaling * Sur_female;		
                       females = females + 1;
                     }	
                   else if (ind.sex =="M"){
                     ind.fitnessScaling = ind.fitnessScaling * Sur_male;	
                     males = males + 1;
                   }
                   
                   maleRatio = males / (males + females);
                   femaleRatio = 1 - maleRatio;
                   
                   sexComp = (Sur_male * maleRatio + Sur_female * femaleRatio);
                   
                   subpop.fitnessScaling = min(K[o] / min(subpop.individualCount *
                                                            mean(age_survival)* sexComp) , Low_density_benefit);	
                   
                   
                 }
                 o = o + 1;
               }
               
             }
  ),
  
  slim_block(modifyChild(), {
    
    if (!child.genome1.containsMarkerMutation(m11, L-1)){
      FMp1 = parent1.genome1.mutationsOfType(m7);
      FMp2 = parent1.genome1.mutationsOfType(m8);
      FMp3 = parent1.genome1.mutationsOfType(m9);
      FMp4 = parent1.genome1.mutationsOfType(m10);
      FMMark = parent1.genome1.mutationsOfType(m11);
      FMMut = parent1.genome1.mutationsOfType(m2);	
      
      child.genome1.addMutations(c(FMp1, FMp2, FMp3,FMp4, FMMark, FMMut));
    }
    
    CMp1 = child.genome2.mutationsOfType(m7);
    CMp2 = child.genome2.mutationsOfType(m8);
    CMp3 = child.genome2.mutationsOfType(m9);
    CMp4 = child.genome2.mutationsOfType(m10);
    CMMark = child.genome2.mutationsOfType(m11);
    CMMut = child.genome2.mutationsOfType(m2);
    
    child.genome2.removeMutations(c(CMp1, CMp2, CMp3, CMp4, CMMark, CMMut));
    
    return(T);
  }
  ),
  
  slim_block(18000, late(), { 
    
    ##nuclear mutations
    nuc  = sim.subpopulations.individuals.genomes.mutationsOfType(m1);
    nuc1 = sim.subpopulations.individuals.genomes.mutationsOfType(m3);
    nuc2 = sim.subpopulations.individuals.genomes.mutationsOfType(m4);
    nuc3 = sim.subpopulations.individuals.genomes.mutationsOfType(m5);
    nuc4 = sim.subpopulations.individuals.genomes.mutationsOfType(m6);
    
    defineGlobal("Nucs", c(nuc, nuc1, nuc2, nuc3, nuc4));	
    
    mit  = sim.subpopulations.individuals.genome1.mutationsOfType(m2); 
    mit1 = sim.subpopulations.individuals.genome1.mutationsOfType(m7);
    mit2 = sim.subpopulations.individuals.genome1.mutationsOfType(m8);
    mit3 = sim.subpopulations.individuals.genome1.mutationsOfType(m9);
    mit4 = sim.subpopulations.individuals.genome1.mutationsOfType(m10);
    
    defineGlobal("Mits", c(mit, mit1, mit2, mit3, mit4));	
    
    Pop1 = p1.sampleIndividuals(30);
    Pop2 = p2.sampleIndividuals(30);
    Pop3 = p3.sampleIndividuals(30);
    Pop4 = p4.sampleIndividuals(30);
    
    ## output FST data with slimr verbs.
    r_output(calcFST(Pop1.genome1,Pop2.genome1, Mits), "FMp12");
    r_output(calcFST(Pop1.genome1,Pop3.genome1, Mits), "FMp13");
    r_output(calcFST(Pop1.genome1,Pop4.genome1, Mits), "FMp14");
    r_output(calcFST(Pop2.genome1,Pop3.genome1, Mits), "FMp23");
    r_output(calcFST(Pop2.genome1,Pop4.genome1, Mits), "FMp24");
    r_output(calcFST(Pop3.genome1,Pop4.genome1, Mits), "FMp34");
    r_output(calcFST(Pop1.genomes,Pop2.genomes, Nucs), "FNp12");
    r_output(calcFST(Pop1.genomes,Pop3.genomes, Nucs), "FNp13");
    r_output(calcFST(Pop1.genomes,Pop4.genomes, Nucs), "FNp14");
    r_output(calcFST(Pop2.genomes,Pop3.genomes, Nucs), "FNp23");
    r_output(calcFST(Pop2.genomes,Pop4.genomes, Nucs), "FNp24");
    r_output(calcFST(Pop3.genomes,Pop4.genomes, Nucs), "FNp34");
    
  } )) -> humpy_script_0

##Slim script for the simulation of the nuclear and mitochondrial DNA. Nm=2
slim_script(
  slim_block(initialize(), 
             {
               ##initializing a nucleotide based non-WF simulation
               initializeSLiMModelType("nonWF");	
               initializeSLiMOptions(nucleotideBased=T);
               
               ##reading a Fasta file as Ancestral string for the individuals
               defineConstant("L", initializeAncestralNucleotides("C:/Users/Domenico/Documents/humpy/Ancestral_Fasta_test.fasta"));
               
               ##Sexual reproduction call
               initializeSex("A");	
               
               ##vector of carrying capacities for each subpopulation and a vector for the benefit of the specific feeding ground of the subpopulation
               defineConstant("K", c(1000.0, 1000.0, 1000.0, 1000.0));	
               
               ##number of populations
               defineConstant("P", 4);
               
               ##Number of individuals per subpopulation at the start of the simulation 
               defineConstant("N", 10);	
               
               ##Age of maturity and maximum age of the population are defined here (will be changed to a sex-based age of maturity later on. This would require additionial changed to the mortality age table
               defineConstant("Age_mat", 5);
               defineConstant("Age_max", 80);
               
               ##to counteract the early boost in fitness from the carrying capacity formula a low-density value is defined so that at low density mortality is only due to age and sex.
               defineConstant("Low_density_benefit", asFloat(1.0));
               
               #Sex based survivability
               defineConstant("Sur_male", 1.0);
               defineConstant("Sur_female", 1.0);
               
               ##constants for migrations
               defineConstant("Nm", 2);
               defineConstant("Gen", 12);	
               
               Mort = float(Age_max+1);
               Mort[(Age_max-10):Age_max] = c(0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);
               Mort[0:Age_mat] = c(0.05, 0.05, 0.025, 0.025, 0.025, 0.025);
               
               defineConstant ("Mort_age", Mort);
               
               for (i in 1:10) initializeMutationTypeNuc(i, 0.5, "f", 0.0);
               initializeMutationType("m11", 0.5, "f", 0.0);
               
               m1.color = "purple";
               m2.color = "pink";
               m3.color = "blue";
               m4.color = "red";
               m5.color = "yellow";
               m6.color = "green";
               m7.color = "blue";
               m8.color = "red";
               m9.color = "yellow";
               m10.color = "green";
               m11.color = "grey";
               
               initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(1.12e-8));	##neutral nuclear DNA region
               initializeGenomicElementType("g2", m2, 1.0, mmJukesCantor(4.3e-6));	##neutral mitochondrial DNA region
               initializeGenomicElement(g1, 0, asInteger(L-22));
               initializeGenomicElement(g2, asInteger(L-21), asInteger(L-1));
               initializeRecombinationRate(c(1e-8, 0), c(asInteger(L-22), asInteger(L-1)));	
               
               g1.color = "blue";
               g2.color = "pink";	
               
             }
  ),
  
  slim_block(reproduction(), 
             
             {	
               if (individual.sex == "F" & individual.age >= Age_mat & individual.tag == 0) {
                 
                 matePop = sample(sim.subpopulations, 1);
                 mate= matePop.sampleIndividuals(1, sex = "M", minAge = Age_mat);
                 
                 ##Checks if the above sample actually contains a male
                 if (size(mate) > 0)
                 { 
                   child = subpop.addCrossed(individual, mate);
                   child.tag = 0;
                 }
                 
                 individual.tag = rdunif(1, min=1, max=2);
               }
               
               ##if the individual has a tag higher than 0 it is lowered by 1 and once the tag reaches zero again the female can reproduce again
               else if (individual.tag > 0) individual.tag = individual.tag - 1;
             }
  ),
  
  slim_block(1, early(),
             {
               for (i in 1:P) sim.addSubpop(i, N);
               
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1N.vcf", m3);
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1M.vcf", m7);
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2N.vcf", m4);	
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2M.vcf", m8);
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3N.vcf", m5);	
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3M.vcf", m9);
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4N.vcf", m6);	
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4M.vcf", m10);
               
               for (subpop in sim.subpopulations) 
               {
                 ind = subpop.individuals;
                 ind.age = rdunif(N, min=0, max=Age_max-10);
                 ind.tag = 0;
               }		
               
               gen1 = sim.subpopulations.individuals;
               gen1.genome1.addNewMutation(m11, 0.0, L-1);
               
             }
  ),
  
  slim_block(early(),
             {
               ##migration event every 20 years
               if (sim.cycle%%Gen ==0) {
                 for (i in 1:Nm)
                 { 
                   migrants = sim.subpopulations.sampleIndividuals(1, sex =  "F");
                   newHomes = c(p2, p3, p4, p1);
                   
                   for (x in 0:(P-1)) 
                   {
                     newHome = newHomes[x];
                     newHome.takeMigrants(migrants[x]);
                   }
                 }
               }
               
               o = 0;
               
               ##life table mortality p1
               for (subpop in sim.subpopulations)
               {
                 if (subpop.individualCount > 0) {
                   inds = subpop.individuals;
                   ages = inds.age;
                   age_mortality = Mort_age[ages];
                   age_survival = 1 - age_mortality;
                   
                   ##fitness is calculated for all individuals based on their age and feeding ground
                   inds.fitnessScaling = age_survival;
                   
                   ##the fitnessScaling is altered by a survival factor based on the sex of the individual
                   females = 0;
                   males = 0;		
                   
                   for (ind in inds)
                     if (ind.sex =="F"){
                       ind.fitnessScaling = ind.fitnessScaling * Sur_female;		
                       females = females + 1;
                     }	
                   else if (ind.sex =="M"){
                     ind.fitnessScaling = ind.fitnessScaling * Sur_male;	
                     males = males + 1;
                   }
                   
                   maleRatio = males / (males + females);
                   femaleRatio = 1 - maleRatio;
                   
                   sexComp = (Sur_male * maleRatio + Sur_female * femaleRatio);
                   
                   subpop.fitnessScaling = min(K[o] / min(subpop.individualCount *
                                                            mean(age_survival) * sexComp) , Low_density_benefit);	
                   
                   
                 }
                 o = o + 1;
               }
               
             }
  ),
  
  slim_block(modifyChild(), {
    
    if (!child.genome1.containsMarkerMutation(m11, L-1)){
      FMp1 = parent1.genome1.mutationsOfType(m7);
      FMp2 = parent1.genome1.mutationsOfType(m8);
      FMp3 = parent1.genome1.mutationsOfType(m9);
      FMp4 = parent1.genome1.mutationsOfType(m10);
      FMMark = parent1.genome1.mutationsOfType(m11);
      FMMut = parent1.genome1.mutationsOfType(m2);	
      
      child.genome1.addMutations(c(FMp1, FMp2, FMp3,FMp4, FMMark, FMMut));
    }
    
    CMp1 = child.genome2.mutationsOfType(m7);
    CMp2 = child.genome2.mutationsOfType(m8);
    CMp3 = child.genome2.mutationsOfType(m9);
    CMp4 = child.genome2.mutationsOfType(m10);
    CMMark = child.genome2.mutationsOfType(m11);
    CMMut = child.genome2.mutationsOfType(m2);
    
    child.genome2.removeMutations(c(CMp1, CMp2, CMp3, CMp4, CMMark, CMMut));
    
    return(T);
  }
  ),
  
  slim_block(18000, late(), { 
    
    ##nuclear mutations
    nuc  = sim.subpopulations.individuals.genomes.mutationsOfType(m1);
    nuc1 = sim.subpopulations.individuals.genomes.mutationsOfType(m3);
    nuc2 = sim.subpopulations.individuals.genomes.mutationsOfType(m4);
    nuc3 = sim.subpopulations.individuals.genomes.mutationsOfType(m5);
    nuc4 = sim.subpopulations.individuals.genomes.mutationsOfType(m6);
    
    defineGlobal("Nucs", c(nuc, nuc1, nuc2, nuc3, nuc4));	
    
    mit  = sim.subpopulations.individuals.genome1.mutationsOfType(m2); 
    mit1 = sim.subpopulations.individuals.genome1.mutationsOfType(m7);
    mit2 = sim.subpopulations.individuals.genome1.mutationsOfType(m8);
    mit3 = sim.subpopulations.individuals.genome1.mutationsOfType(m9);
    mit4 = sim.subpopulations.individuals.genome1.mutationsOfType(m10);
    
    defineGlobal("Mits", c(mit, mit1, mit2, mit3, mit4));	
    
    Pop1 = p1.sampleIndividuals(30);
    Pop2 = p2.sampleIndividuals(30);
    Pop3 = p3.sampleIndividuals(30);
    Pop4 = p4.sampleIndividuals(30);
    
    ## output FST data with slimr verbs.
    r_output(calcFST(Pop1.genome1,Pop2.genome1, Mits), "FMp12");
    r_output(calcFST(Pop1.genome1,Pop3.genome1, Mits), "FMp13");
    r_output(calcFST(Pop1.genome1,Pop4.genome1, Mits), "FMp14");
    r_output(calcFST(Pop2.genome1,Pop3.genome1, Mits), "FMp23");
    r_output(calcFST(Pop2.genome1,Pop4.genome1, Mits), "FMp24");
    r_output(calcFST(Pop3.genome1,Pop4.genome1, Mits), "FMp34");
    r_output(calcFST(Pop1.genomes,Pop2.genomes, Nucs), "FNp12");
    r_output(calcFST(Pop1.genomes,Pop3.genomes, Nucs), "FNp13");
    r_output(calcFST(Pop1.genomes,Pop4.genomes, Nucs), "FNp14");
    r_output(calcFST(Pop2.genomes,Pop3.genomes, Nucs), "FNp23");
    r_output(calcFST(Pop2.genomes,Pop4.genomes, Nucs), "FNp24");
    r_output(calcFST(Pop3.genomes,Pop4.genomes, Nucs), "FNp34");
    
  } )) -> humpy_script_2

##Slim script for the simulation of the nuclear and mitochondrial DNA. Nm=10
slim_script(
  slim_block(initialize(), 
             {
               ##initializing a nucleotide based non-WF simulation
               initializeSLiMModelType("nonWF");	
               initializeSLiMOptions(nucleotideBased=T);
               
               ##reading a Fasta file as Ancestral string for the individuals
               defineConstant("L", initializeAncestralNucleotides("C:/Users/Domenico/Documents/humpy/Ancestral_Fasta_test.fasta"));
               
               ##Sexual reproduction call
               initializeSex("A");	
               
               ##vector of carrying capacities for each subpopulation and a vector for the benefit of the specific feeding ground of the subpopulation
               defineConstant("K", c(1000.0, 1000.0, 1000.0, 1000.0));
               
               ##number of populations
               defineConstant("P", 4);
               
               ##Number of individuals per subpopulation at the start of the simulation 
               defineConstant("N", 10);	
               
               ##Age of maturity and maximum age of the population are defined here (will be changed to a sex-based age of maturity later on. This would require additionial changed to the mortality age table
               defineConstant("Age_mat", 5);
               defineConstant("Age_max", 80);
               
               ##to counteract the early boost in fitness from the carrying capacity formula a low-density value is defined so that at low density mortality is only due to age and sex.
               defineConstant("Low_density_benefit", asFloat(1.0));
               
               #Sex based survivability
               defineConstant("Sur_male", 1.0);
               defineConstant("Sur_female", 1.0);
               
               ##constants for migrations
               defineConstant("Nm", 10);
               defineConstant("Gen", 12);	
               
               Mort = float(Age_max+1);
               Mort[(Age_max-10):Age_max] = c(0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);
               Mort[0:Age_mat] = c(0.05, 0.05, 0.025, 0.025, 0.025, 0.025);
               
               defineConstant ("Mort_age", Mort);
               
               for (i in 1:10) initializeMutationTypeNuc(i, 0.5, "f", 0.0);
               initializeMutationType("m11", 0.5, "f", 0.0);
               
               m1.color = "purple";
               m2.color = "pink";
               m3.color = "blue";
               m4.color = "red";
               m5.color = "yellow";
               m6.color = "green";
               m7.color = "blue";
               m8.color = "red";
               m9.color = "yellow";
               m10.color = "green";
               m11.color = "grey";
               
               initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(1.12e-8));	##neutral nuclear DNA region
               initializeGenomicElementType("g2", m2, 1.0, mmJukesCantor(4.3e-6));	##neutral mitochondrial DNA region 
               initializeGenomicElement(g1, 0, asInteger(L-22));
               initializeGenomicElement(g2, asInteger(L-21), asInteger(L-1));
               initializeRecombinationRate(c(1e-8, 0), c(asInteger(L-22), asInteger(L-1)));	
               
               g1.color = "blue";
               g2.color = "pink";	
               
             }
  ),
  
  slim_block(reproduction(), 
             
             {	
               if (individual.sex == "F" & individual.age >= Age_mat & individual.tag == 0) {
                 
                 matePop = sample(sim.subpopulations, 1);
                 mate= matePop.sampleIndividuals(1, sex = "M", minAge = Age_mat);
                 
                 ##Checks if the above sample actually contains a male
                 if (size(mate) > 0)
                 { 
                   child = subpop.addCrossed(individual, mate);
                   child.tag = 0;
                 }
                 
                 individual.tag = rdunif(1, min=1, max=2);
               }
               
               ##if the individual has a tag higher than 0 it is lowered by 1 and once the tag reaches zero again the female can reproduce again
               else if (individual.tag > 0) individual.tag = individual.tag - 1;
             }
  ),
  
  slim_block(1, early(),
             {
               for (i in 1:P) sim.addSubpop(i, N);
               
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1N.vcf", m3);
               p1.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group1M.vcf", m7);
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2N.vcf", m4);	
               p2.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group2M.vcf", m8);
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3N.vcf", m5);	
               p3.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group3M.vcf", m9);
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4N.vcf", m6);	
               p4.genomes.readFromVCF("C:/Users/Domenico/Documents/humpy/group4M.vcf", m10);
               
               for (subpop in sim.subpopulations) 
               {
                 ind = subpop.individuals;
                 ind.age = rdunif(N, min=0, max=Age_max-10);
                 ind.tag = 0;
               }		
               
               gen1 = sim.subpopulations.individuals;
               gen1.genome1.addNewMutation(m11, 0.0, L-1);
               
             }
  ),
  
  slim_block(early(),
             {
               ##migration event every 20 years
               if (sim.cycle%%Gen ==0) {
                 for (i in 1:Nm)
                 { 
                   migrants = sim.subpopulations.sampleIndividuals(1, sex =  "F");
                   newHomes = c(p2, p3, p4, p1);
                   
                   for (x in 0:(P-1))
                   {
                     newHome = newHomes[x];
                     newHome.takeMigrants(migrants[x]);
                   }
                 }
               }
               
               o = 0;
               
               ##life table mortality p1
               for (subpop in sim.subpopulations)
               {
                 if (subpop.individualCount > 0) {
                   inds = subpop.individuals;
                   ages = inds.age;
                   age_mortality = Mort_age[ages];
                   age_survival = 1 - age_mortality;
                   
                   ##fitness is calculated for all individuals based on their age and feeding ground
                   inds.fitnessScaling = age_survival;
                   
                   ##the fitnessScaling is altered by a survival factor based on the sex of the individual
                   females = 0;
                   males = 0;		
                   
                   for (ind in inds)
                     if (ind.sex =="F"){
                       ind.fitnessScaling = ind.fitnessScaling * Sur_female;		
                       females = females + 1;
                     }	
                   else if (ind.sex =="M"){
                     ind.fitnessScaling = ind.fitnessScaling * Sur_male;	
                     males = males + 1;
                   }
                   
                   maleRatio = males / (males + females);
                   femaleRatio = 1 - maleRatio;
                   
                   sexComp = (Sur_male * maleRatio + Sur_female * femaleRatio);
                   
                   subpop.fitnessScaling = min(K[o] / min(subpop.individualCount *
                                                            mean(age_survival) * sexComp) , Low_density_benefit);	
                   
                   
                 }
                 o = o + 1;
               }
               
             }
  ),
  
  slim_block(modifyChild(), {
    
    if (!child.genome1.containsMarkerMutation(m11, L-1)){
      FMp1 = parent1.genome1.mutationsOfType(m7);
      FMp2 = parent1.genome1.mutationsOfType(m8);
      FMp3 = parent1.genome1.mutationsOfType(m9);
      FMp4 = parent1.genome1.mutationsOfType(m10);
      FMMark = parent1.genome1.mutationsOfType(m11);
      FMMut = parent1.genome1.mutationsOfType(m2);	
      
      child.genome1.addMutations(c(FMp1, FMp2, FMp3,FMp4, FMMark, FMMut));
    }
    
    CMp1 = child.genome2.mutationsOfType(m7);
    CMp2 = child.genome2.mutationsOfType(m8);
    CMp3 = child.genome2.mutationsOfType(m9);
    CMp4 = child.genome2.mutationsOfType(m10);
    CMMark = child.genome2.mutationsOfType(m11);
    CMMut = child.genome2.mutationsOfType(m2);
    
    child.genome2.removeMutations(c(CMp1, CMp2, CMp3, CMp4, CMMark, CMMut));
    
    return(T);
  }
  ),
  
  slim_block(18000, late(), { 
    
    ##nuclear mutations
    nuc  = sim.subpopulations.individuals.genomes.mutationsOfType(m1);
    nuc1 = sim.subpopulations.individuals.genomes.mutationsOfType(m3);
    nuc2 = sim.subpopulations.individuals.genomes.mutationsOfType(m4);
    nuc3 = sim.subpopulations.individuals.genomes.mutationsOfType(m5);
    nuc4 = sim.subpopulations.individuals.genomes.mutationsOfType(m6);
    
    defineGlobal("Nucs", c(nuc, nuc1, nuc2, nuc3, nuc4));	
    
    mit  = sim.subpopulations.individuals.genome1.mutationsOfType(m2); 
    mit1 = sim.subpopulations.individuals.genome1.mutationsOfType(m7);
    mit2 = sim.subpopulations.individuals.genome1.mutationsOfType(m8);
    mit3 = sim.subpopulations.individuals.genome1.mutationsOfType(m9);
    mit4 = sim.subpopulations.individuals.genome1.mutationsOfType(m10);
    
    defineGlobal("Mits", c(mit, mit1, mit2, mit3, mit4));	
    
    Pop1 = p1.sampleIndividuals(30);
    Pop2 = p2.sampleIndividuals(30);
    Pop3 = p3.sampleIndividuals(30);
    Pop4 = p4.sampleIndividuals(30);
    
    ## output FST data with slimr verbs.
    r_output(calcFST(Pop1.genome1,Pop2.genome1, Mits), "FMp12");
    r_output(calcFST(Pop1.genome1,Pop3.genome1, Mits), "FMp13");
    r_output(calcFST(Pop1.genome1,Pop4.genome1, Mits), "FMp14");
    r_output(calcFST(Pop2.genome1,Pop3.genome1, Mits), "FMp23");
    r_output(calcFST(Pop2.genome1,Pop4.genome1, Mits), "FMp24");
    r_output(calcFST(Pop3.genome1,Pop4.genome1, Mits), "FMp34");
    r_output(calcFST(Pop1.genomes,Pop2.genomes, Nucs), "FNp12");
    r_output(calcFST(Pop1.genomes,Pop3.genomes, Nucs), "FNp13");
    r_output(calcFST(Pop1.genomes,Pop4.genomes, Nucs), "FNp14");
    r_output(calcFST(Pop2.genomes,Pop3.genomes, Nucs), "FNp23");
    r_output(calcFST(Pop2.genomes,Pop4.genomes, Nucs), "FNp24");
    r_output(calcFST(Pop3.genomes,Pop4.genomes, Nucs), "FNp34");
    
  } )) -> humpy_script_10