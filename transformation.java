import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;



public class transformation {
    /**
     * @param args
     */
    public static void main(String[] args) {
       
        double [] param=new double[16] ;
        param[0]=Integer.parseInt(args[0],10);
        param[1]=Double.parseDouble(args[1]);
        param[2]=Double.parseDouble(args[2]);
        param[3]=Double.parseDouble(args[3]);
        param[4]=Double.parseDouble(args[4]);
        param[5]=Double.parseDouble(args[5]);
        param[6]=Double.parseDouble(args[6]);
        param[7]=Double.parseDouble(args[7]);
        param[8]=Double.parseDouble(args[8]);
        param[9]=Double.parseDouble(args[9]);
        param[10]=Double.parseDouble(args[10]);
        param[11]=Integer.parseInt(args[11],10);
        param[12]=Integer.parseInt(args[12],10);
        param[13]=Integer.parseInt(args[13],10);
        param[14]=Integer.parseInt(args[14],10);
        param[15]=Double.parseDouble(args[15]);
        String rectemp="recombinant"+Integer.parseInt(args[0],10);
        PrintWriter recombinant = MakeFile(rectemp);
       
        String nuttemp="nutrient"+Integer.parseInt(args[0],10);
        PrintWriter nutrient = MakeFile(nuttemp);
       
        //String inctemp="/Users/danesh/Documents/workspace/transformation/noncompetent.txt";
        String inctemp="noncompetent"+Integer.parseInt(args[0],10);
        PrintWriter noncompetent = MakeFile(inctemp);
       
        //String fragtemp="/Users/danesh/Documents/workspace/transformation/fragment.txt";
        String fragtemp="fragment"+Integer.parseInt(args[0],10);
        PrintWriter genfrag = MakeFile(fragtemp);
       
       
       
       
       
       
        int replicate=1;
       
       
        for(int it=0;it<replicate;it++){
       
       
       
       
        int loci=Integer.parseInt(args[14],10);
        double DNAsize= param[1];
        int dimen=(int) Math.pow(2,loci);
        int dimens=3* (int) Math.pow(2,loci);
        int dimenfr=(int) Math.pow(2,DNAsize);
        int DNAdimen=(int) Math.pow(2,DNAsize)*(loci-(int) DNAsize+1);
   
        double taum=0.003;
        int cut=10;
        double time=0;
        Date d = new Date();
        RandomEngine engine = new DRand((int) d.getTime());
       
       
        double recrate= param[15];//0.000001;
        double uptakerate=param[2];//0.0000001;//;0.000001;//0.0000001;
        double foodvalue=param[3];//0.000001;//0.00001;//0.1-> 0.01
        double degredation=param[4];//1;//10;
        double upper=10000;
        double selectionco=param[5];//0.02;//0.01;
        double mutation=param[6];//0.000000001;//0.000001;
        double populationsize=0;
        double fragmentsize=0;
        double capacity=param[7];//100000000;//9757186
        double growthrate=param[8];//0.1;
        double competentsize=0;
        double cost= param[9];//0.0000001;//0.00001;//0.00001;
        double initialsize=0;
        double deathrate=param[10];// 0.0001;
       
        int [] sizen= new int[(int) Math.pow(2,loci)];
        int [] sizer= new int[(int) Math.pow(2,loci)];
        int [] sizeu= new int[(int) Math.pow(2,loci)];
       
        for(int i=0;i<dimen;i++){sizen[i]=0;sizer[i]=0;sizeu[i]=0;}
    //    Integer.parseInt(args[11],10);
    //    param[12]=Integer.parseInt(args[12],10);
    //    param[13]=Integer.parseInt(args[13],10);
        sizen[0]=Integer.parseInt(args[11],10);//10000;
        sizer[0]=Integer.parseInt(args[12],10);
        sizeu[0]=Integer.parseInt(args[13],10);
        for(int i=0;i<dimen;i++){populationsize=populationsize+sizen[i]+sizer[i]+sizeu[i];}
        for(int i=0;i<dimen;i++){initialsize=initialsize+sizen[i]+sizer[i]+sizeu[i];}
        for(int i=0;i<dimen;i++){competentsize=competentsize+sizen[i]+sizer[i];}
        double [] fragment=new double[DNAdimen];
        for(int r=0;r<fragment.length;r++){fragment[r]=0;}
       
    //    System.out.println(Arrays.toString(sizen));
    //    System.out.println(Arrays.toString(sizer));
    //    System.out.println(Arrays.toString(sizeu));
        //////////////////////////////////////////////birth
        double[] birthn= new double[dimen];
        double[] birthr= new double[dimen];
        double[] birthu= new double[dimen];
        for(int i=0;i<dimen;i++){
             birthn[i]=sizen[i]*(growthrate+uptakerate*foodvalue*fragmentsize-cost)*(1-(populationsize/capacity));
            birthr[i]=sizer[i]*(growthrate+uptakerate*foodvalue*fragmentsize-cost)*(1-(populationsize/capacity));   
            birthu[i]=sizeu[i]*(growthrate)*(1-(populationsize/capacity));   
        }

        double [] Fitness= new double[(int) Math.pow(2,loci)];
        Fitness=Fitness(selectionco,loci,dimen);
        //System.out.println(Arrays.toString(Fitness));
       
       
        //Definition matrix death rate
        double [][] deathratem=new double[dimen][4*(loci-(int) DNAsize+1)];
        for(int y=0;y<dimen;y++){
            String holder=ToBinary(y,loci);
            for(int u=0;u<4*(loci-(int) DNAsize+1);u=u+4){
                String subholder=holder.substring(u/4, u/4+ (int) DNAsize);
                deathratem[y][u]=Integer.parseInt(subholder, 2);
                deathratem[y][u+1]=sizen[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
                deathratem[y][u+2]=sizer[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
                deathratem[y][u+3]=sizeu[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
            }
        }

       
        ///////////////////////////////////////////////////////////////////mutation
        ArrayList <double []> eventmutn=new ArrayList <double []> ();
        ArrayList <double []> eventmutr=new ArrayList <double []> ();
        ArrayList <double []> eventmutu=new ArrayList <double []> ();
        for(int i=0;i<dimen;i++){
            for(int j=0;j<dimen;j++){
                hamming h= new hamming (ToBinary(i,loci),ToBinary(j,loci));
                if(h.getHammingDistance()==1){
                    double [] temp={Double.parseDouble(Integer.toString(i)),Double.parseDouble(Integer.toString(j)),mutation*sizen[i]};
                    eventmutn.add(temp);
                    double [] temp1={Double.parseDouble(Integer.toString(i)),Double.parseDouble(Integer.toString(j)),mutation*sizer[i]};
                    eventmutr.add(temp1);
                    double [] temp2={Double.parseDouble(Integer.toString(i)),Double.parseDouble(Integer.toString(j)),mutation*sizeu[i]};
                    eventmutu.add(temp2);
                    };   
            };
        };
        double [][] mutn=new double[eventmutn.size()][];
        double [][] mutr=new double[eventmutn.size()][];
        double [][] mutu=new double[eventmutn.size()][];
        eventmutn.toArray(mutn);
        eventmutr.toArray(mutr);
        eventmutu.toArray(mutu);
        ///////////////////////////////////////////////////////////////////
       
        //printmatrix(mutn);
       
       
       
        //////////////////////////////////////////////////////////degredation
        double [] degred=new double[DNAdimen];
        for(int r=0;r<fragment.length;r++){
            degred[r]=degredation*fragment[r];
            }
        ///////////////////////////////////////////////////////////
       
        //////////////////////////////////////////////////////////uptake
        double [][] digest=new double[dimen][DNAdimen];
        for(int i=0;i<dimen;i++){
            for(int j=0;j<(int) DNAdimen;j++){
                digest[i][j]=sizen[i]*uptakerate*fragment[j];       
            }
        }
        /////////////////////////////////////////////////////////
       
       
       
        //////////////////////////////////////////////////////////uptake+recombination
        double [][][] recombinationrate=new double[dimen][DNAdimen][2];       
        for(int i=0;i<dimen;i++){
            String holder=ToBinary(i,loci);
            int counter=0;
            for(int j=0;j<(int) DNAdimen;j++){
                String tempor=ToBinary(counter,(int) DNAsize);
                int start=j/dimenfr;
                String third="";
                third=third.concat(holder.substring(0, start).concat(tempor).concat(holder.substring(start+(int) DNAsize, loci)));
                 recombinationrate[i][j][0]=Integer.parseInt(third, 2);
                 recombinationrate[i][j][1]=sizer[i]*uptakerate*fragment[j];       
                counter++;
                if(counter%dimenfr==0){counter=0;}
            }
        }
        /////////////////////////////////////////////////////////
        //Start of the loop
        ///////////////////////////////////////////////////////
       
        int counterloop=0;
       
        int entry=0;
       
        while(true){
            double tau=taum;
           
           
            if(counterloop%5000==0){
                entry++;
               
                recombinant.println(time);
                recombinant.println(Arrays.toString(sizer));
               
                nutrient.println(time);
                nutrient.println(Arrays.toString(sizen));
               
                noncompetent.println(time);
                noncompetent.println(Arrays.toString(sizeu));
               
                genfrag.println(time);
                genfrag.println(Arrays.toString(fragment));
               

            }
            counterloop++;
        ArrayList <double []> majev=new ArrayList <double []> ();
        ArrayList <int []> majid=new ArrayList <int []> ();
       
       
        ArrayList <double []> minev=new ArrayList <double []> ();
        ArrayList <int []> minid=new ArrayList <int []> ();
       
       
        for(int j=0;j<dimen;j++){
            if(sizen[j]>cut){int [] temp={1,j};majid.add(temp);} else {int [] temp={1,j};minid.add(temp);};
            if(sizer[j]>cut){int [] temp={2,j};majid.add(temp);} else {int [] temp={2,j};minid.add(temp);};
            if(sizeu[j]>cut){int [] temp={3,j};majid.add(temp);} else {int [] temp={3,j};minid.add(temp);};
        }
       
        for(int j=0;j<fragment.length;j++){
            if(fragment[j]>=cut){int [] temp={4,j};majid.add(temp);} else {int [] temp={4,j};minid.add(temp);};
        }
       
       

       
        for(int i=0;i<dimen;i++){
            if(checking( majid,new int[] { 1, i})){majev.add(new double[] {1,1,i,birthn[i],1,1,1});} else if(checking( minid,new int[] { 1, i})) {minev.add(new double[] {1,1,i,birthn[i],1,1,1});};
            if(checking( majid,new int[] { 2, i})){majev.add(new double[] {1,2,i,birthr[i],1,1,1});} else if(checking( minid,new int[] { 2, i})) {minev.add(new double[] {1,2,i,birthr[i],1,1,1});};
            if(checking( majid,new int[] { 3, i})){majev.add(new double[] {1,3,i,birthu[i],1,1,1});} else if(checking( minid,new int[] { 3, i})) {minev.add(new double[] {1,3,i,birthu[i],1,1,1});};   
        }
       
       
       
           
        for(int y=0;y<dimen;y++){
            for(int u=0;u<4*(loci-(int) DNAsize+1);u=u+4){
                if(checking( majid,new int[] { 1, y})&&checking( majid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]})){
                    majev.add(new double[] {2,1,y,deathratem[y][u+1],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+1]});
                }
                else
                    if(checking( minid,new int[] { 1, y})||checking( minid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]}))
                {
                    minev.add(new double[] {2,1,y,deathratem[y][u+1],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+1]});   
                }
               
               
                if(checking( majid,new int[] { 2, y})&&checking( majid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]})){
                    majev.add(new double[] {2,2,y,deathratem[y][u+2],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+2]});
                }
                else
                    if(checking( minid,new int[] { 2, y})||checking( minid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]}))
                {
                    minev.add(new double[] {2,2,y,deathratem[y][u+2],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+2]});   
                }
                   
                if(checking( majid,new int[] { 3, y})&&checking( majid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]})){
                    majev.add(new double[] {2,3,y,deathratem[y][u+3],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+3]});
                }
                else
                    if(checking( minid,new int[] { 3, y})||checking( minid,new int[] {4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u]})){
                    minev.add(new double[] {2,3,y,deathratem[y][u+3],4,(int) Math.pow(2,DNAsize)*(u/4)+(int) deathratem[y][u],deathratem[y][u+3]});   
                }
       
            }
        }
       
        for(int i=0;i<mutn.length;i++){
            if(checking( majid,new int[] { 1, (int) mutn[i][0]})&&checking( majid,new int[] {1,(int) mutn[i][1]})){
                majev.add(new double[] {3,1,mutn[i][0],mutn[i][2],1,mutn[i][1],mutn[i][2]});
            }
            else
                if(checking( minid,new int[] { 1, (int) mutn[i][0]})||checking( minid,new int[] {1,(int) mutn[i][1]}))
            {
                minev.add(new double[] {3,1,mutn[i][0],mutn[i][2],1,mutn[i][1],mutn[i][2]});   
            }
           
           
            if(checking( majid,new int[] { 2, (int) mutr[i][0]})&&checking( majid,new int[] {2,(int) mutr[i][1]})){
                majev.add(new double[] {3,2,mutr[i][0],mutr[i][2],2,mutr[i][1],mutr[i][2]});
            }
            else
                if(checking( minid,new int[] { 2, (int) mutr[i][0]})||checking( minid,new int[] {2,(int) mutr[i][1]}))
            {
                minev.add(new double[] {3,2,mutr[i][0],mutr[i][2],2,mutr[i][1],mutr[i][2]});   
            }
           
           
            if(checking( majid,new int[] { 3, (int) mutu[i][0]})&& checking( majid,new int[] {3,(int) mutu[i][1]})){
                majev.add(new double[] {3,3,mutu[i][0],mutu[i][2],3,mutu[i][1],mutu[i][2]});
            }
            else
                if(checking( minid,new int[] { 3, (int) mutu[i][0]})||checking( minid,new int[] {3,(int) mutu[i][1]}))
            {
                minev.add(new double[] {3,3,mutu[i][0],mutu[i][2],3,mutu[i][1],mutu[i][2]});   
            }
        }
       
       
        //eshtebah dar nazar begir kelasi ke tolid mishavad
        for(int i=0;i<dimen;i++){
            for(int j=0;j<(int) DNAdimen;j++){
                if(checking( majid,new int[] { 1, i})&&checking( majid,new int[] {4,j})){
                    majev.add(new double[] {4,1,i,digest[i][j],4,j,digest[i][j]});
                }
                else
                    if(checking( minid,new int[] { 1, i})||checking( minid,new int[] {4,j})){
                    minev.add(new double[] {4,1,i,digest[i][j],4,j,digest[i][j]});   
                }
               
               
                if(checking( majid,new int[] { 2, i})&& checking( majid,new int[] {4,j})&& checking( majid,new int[] {2,(int) recombinationrate[i][j][0]})){
                    majev.add(new double[] {5,2,i,recombinationrate[i][j][1],recombinationrate[i][j][0],j,recombinationrate[i][j][1]});
                }
                else
                    if(checking( minid,new int[] { 2, i})||checking( minid,new int[] {4,j})||checking( minid,new int[] {2,(int) recombinationrate[i][j][0]})){
                    minev.add(new double[] {5,2,i,recombinationrate[i][j][1],recombinationrate[i][j][0],j,recombinationrate[i][j][1]});   
                }   
            }
        }
       
       
        for(int i=0;i<degred.length;i++){
            if(checking( majid,new int[] {4, i})){majev.add(new double[] {6,4,i,degred[i],1,1,1});}
            else if(checking( minid,new int[] {4, i})){minev.add(new double[] {6,4,i,degred[i],1,1,1});}
        }
       
       
        double sumevents=0;
        double flag=0;
        for(int i=0;i<minev.size();i++){sumevents=sumevents+minev.get(i)[3];if(minev.get(i)[3]<0){flag=minev.get(i)[0];System.out.println(minev.get(i)[1]);System.out.println(minev.get(i)[2]);System.out.println(minev.get(i)[3]);System.out.println(Arrays.toString(birthn));break;}};
        if(flag!=0){System.out.println(flag);break;}
        Exponential expon = new Exponential(sumevents, engine);
        double nextevent= expon.nextDouble();
       
       
   
       

        if(nextevent<taum&&nextevent>0){
       
            //minor events
       
            double tem=0;
            double pemp=sumevents*engine.nextDouble();
            int counter=0;
            try{
                do{
                    tem=minev.get(counter)[3]+tem;
                    counter++;
                }while(tem<=pemp);
                counter--;
            }
            catch(Exception e){System.out.println(e);}
       
       
       
            if(minev.get(counter)[0]==1){
                if(minev.get(counter)[1]==1){sizen[(int) minev.get(counter)[2]]=sizen[(int) minev.get(counter)[2]]+1;};
                if(minev.get(counter)[1]==2){sizer[(int) minev.get(counter)[2]]=sizer[(int) minev.get(counter)[2]]+1;};
                if(minev.get(counter)[1]==3){sizeu[(int) minev.get(counter)[2]]=sizeu[(int) minev.get(counter)[2]]+1;};
            };
       
       
            if(minev.get(counter)[0]==2){
                if(minev.get(counter)[1]==1){
                    sizen[(int) minev.get(counter)[2]]=sizen[(int) minev.get(counter)[2]]-1;
                    fragment[(int) minev.get(counter)[5]]=fragment[(int) minev.get(counter)[5]]+1;
                };
           
                if(minev.get(counter)[1]==2){
                    sizer[(int) minev.get(counter)[2]]=sizer[(int) minev.get(counter)[2]]-1;
                    fragment[(int) minev.get(counter)[5]]=fragment[(int) minev.get(counter)[5]]+1;
                };
           
                if(minev.get(counter)[1]==3){
                    sizeu[(int) minev.get(counter)[2]]=sizeu[(int) minev.get(counter)[2]]-1;
                    fragment[(int) minev.get(counter)[5]]=fragment[(int) minev.get(counter)[5]]+1;
                };
            }

            if(minev.get(counter)[0]==3){
                if(minev.get(counter)[1]==1){
                    sizen[(int) minev.get(counter)[2]]=sizen[(int) minev.get(counter)[2]]-1;
                    sizen[(int) minev.get(counter)[5]]=sizen[(int) minev.get(counter)[5]]+1;
                };
           
                if(minev.get(counter)[1]==2){
                    sizer[(int) minev.get(counter)[2]]=sizer[(int) minev.get(counter)[2]]-1;
                    sizer[(int) minev.get(counter)[5]]=sizer[(int) minev.get(counter)[5]]+1;
                };
           
                if(minev.get(counter)[1]==3){
                    sizeu[(int) minev.get(counter)[2]]=sizeu[(int) minev.get(counter)[2]]-1;
                    sizeu[(int) minev.get(counter)[5]]=sizeu[(int) minev.get(counter)[5]]+1;
                };
            }
   
            if(minev.get(counter)[0]==4){
                if(minev.get(counter)[1]==1){
                    fragment[(int) minev.get(counter)[5]]=fragment[(int) minev.get(counter)[5]]-1;
                };
            }
       
            if(minev.get(counter)[0]==5){
                if(minev.get(counter)[1]==2){
                    if(engine.nextDouble()<= recrate){
                    sizer[(int) minev.get(counter)[2]]=sizer[(int) minev.get(counter)[2]]-1;
                    sizer[(int) minev.get(counter)[4]]=sizer[(int) minev.get(counter)[4]]+1;}
                    fragment[(int) minev.get(counter)[5]]=fragment[(int) minev.get(counter)[5]]-1;
                };
            }
       
            if(minev.get(counter)[0]==6){
                if(minev.get(counter)[1]==4){
                    fragment[(int) minev.get(counter)[2]]=fragment[(int) minev.get(counter)[2]]-1;
                };
            }
        }
       
        if(nextevent <= taum){tau=nextevent;}else{tau=taum;}
   
        //major
        double flag2=0;
        for(int i=0;i<majev.size();i++){
            if(majev.get(i)[0]==1){
                if(majev.get(i)[1]==1){if(majev.get(i)[3]>0){sizen[(int) majev.get(i)[2]]=sizen[(int) majev.get(i)[2]]+ new Poisson(majev.get(i)[3]*tau, engine).nextInt();}};
                if(majev.get(i)[1]==2){if(majev.get(i)[3]>0){sizer[(int) majev.get(i)[2]]=sizer[(int) majev.get(i)[2]]+ new Poisson(majev.get(i)[3]*tau, engine).nextInt();}};
                if(majev.get(i)[1]==3){if(majev.get(i)[3]>0){sizeu[(int) majev.get(i)[2]]=sizeu[(int) majev.get(i)[2]]+ new Poisson(majev.get(i)[3]*tau, engine).nextInt();}};
            };
           
           
            if(majev.get(i)[0]==2){
                if(majev.get(i)[1]==1){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizen[(int) majev.get(i)[2]]=sizen[(int) majev.get(i)[2]]-hold;
                    fragment[(int) majev.get(i)[5]]=fragment[(int) majev.get(i)[5]]+hold;}
                };
               
                if(majev.get(i)[1]==2){
                   
                    if(majev.get(i)[3]>0){
                       
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizer[(int) majev.get(i)[2]]=sizer[(int) majev.get(i)[2]]-hold;
                    fragment[(int) majev.get(i)[5]]=fragment[(int) majev.get(i)[5]]+hold;}
                };
               
                if(majev.get(i)[1]==3){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizeu[(int) majev.get(i)[2]]=sizeu[(int) majev.get(i)[2]]-hold;
                    fragment[(int) majev.get(i)[5]]=fragment[(int) majev.get(i)[5]]+hold;}
                };
            }

            if(majev.get(i)[0]==3){
                if(majev.get(i)[1]==1){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizen[(int) majev.get(i)[2]]=sizen[(int) majev.get(i)[2]]-hold;
                    sizen[(int) majev.get(i)[5]]=sizen[(int) majev.get(i)[5]]+hold;}
                };
               
                if(majev.get(i)[1]==2){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizer[(int) majev.get(i)[2]]=sizer[(int) majev.get(i)[2]]-hold;
                    sizer[(int) majev.get(i)[5]]=sizer[(int) majev.get(i)[5]]+hold;}
                };
               
                if(majev.get(i)[1]==3){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    sizeu[(int) majev.get(i)[2]]=sizeu[(int) majev.get(i)[2]]-hold;
                    sizeu[(int) majev.get(i)[5]]=sizeu[(int) majev.get(i)[5]]+hold;}
                };
            }
       
            if(majev.get(i)[0]==4){
                if(majev.get(i)[1]==1){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    if(fragment[(int) majev.get(i)[5]]<hold){flag2=1;System.out.println(hold);System.out.println(fragment[(int) majev.get(i)[5]]);System.out.println(majev.get(i)[2]);System.out.println(majev.get(i)[3]);printmatrix(deathratem);System.out.println(Arrays.toString(birthn));System.out.println("man amadeam");break;}
                    fragment[(int) majev.get(i)[5]]=fragment[(int) majev.get(i)[5]]-hold;}
                };
            }
           
            if(majev.get(i)[0]==5){
                if(majev.get(i)[1]==2){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    if(engine.nextDouble() <= recrate){
                    sizer[(int) majev.get(i)[2]]=sizer[(int) majev.get(i)[2]]-hold;
                    sizer[(int) majev.get(i)[4]]=sizer[(int) majev.get(i)[4]]+hold;}
                    fragment[(int) majev.get(i)[5]]=fragment[(int) majev.get(i)[5]]-hold;}
                };
            }
           
            if(majev.get(i)[0]==6){
                if(majev.get(i)[1]==4){
                    if(majev.get(i)[3]>0){
                    int hold=new Poisson( majev.get(i)[3]*tau, engine).nextInt();
                    fragment[(int) majev.get(i)[2]]=fragment[(int) majev.get(i)[2]]-hold;}
                };
            }
       
       
        }
        if(flag2!=0){break;}
       
        time=time+tau;
        if(time > upper){break;}
       
        populationsize=0;
        for(int i=0;i<dimen;i++){populationsize=populationsize+sizen[i]+sizer[i]+sizeu[i];}
        competentsize=0;
        for(int i=0;i<dimen;i++){competentsize=competentsize+sizen[i]+sizer[i];}
        fragmentsize=0;
        for(int i=0;i<DNAdimen;i++){fragmentsize=fragmentsize+(int) fragment[i];}
       
       
        for(int i=0;i<dimen;i++){
            birthn[i]=sizen[i]*(growthrate+uptakerate*foodvalue*fragmentsize-cost)*(1-populationsize/capacity);
            birthr[i]=sizer[i]*(growthrate+uptakerate*foodvalue*fragmentsize-cost)*(1-populationsize/capacity);   
            birthu[i]=sizeu[i]*(growthrate)*(1-(populationsize/capacity));   
        }
       
       
       
       
        for(int i=0;i<mutn.length;i++){
            mutn[i][2]=sizen[(int) mutn[i][0]]*mutation;
            mutr[i][2]=sizer[(int) mutr[i][0]]*mutation;
            mutu[i][2]=sizeu[(int) mutu[i][0]]*mutation;
        }
       
        for(int y=0;y<dimen;y++){
            for(int u=0;u<4*(loci-(int) DNAsize+1);u=u+4){
                deathratem[y][u+1]=sizen[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
                deathratem[y][u+2]=sizer[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
                deathratem[y][u+3]=sizeu[y]*(deathrate+Fitness[y])/(loci-(int) DNAsize+1);
            }
        }
       
        double fragtot=0;
        for(int i=0;i<(int) DNAdimen;i++){
            fragtot=fragtot+fragment[i];
        }
       
        for(int i=0;i<dimen;i++){
            for(int j=0;j<(int) DNAdimen;j++){
                digest[i][j]=sizen[i]*uptakerate*fragment[j];   
            }
        }
       
        ///eslah degredation rate
        for(int i=0;i<degred.length;i++){
            degred[i]=fragment[i]*degredation;
        }
       
       
        for(int i=0;i<dimen;i++){
            for(int j=0;j<(int) DNAdimen;j++){
                 recombinationrate[i][j][1]=sizer[i]*uptakerate*fragment[j];   
            }
        }
       

    }

        recombinant.println("111111");
        nutrient.println("111111");
        noncompetent.println("111111");
        genfrag.println("111111");

        }

   
    }
   
   
    public static void printmatrix(double [][] arg) {
        for(int i=0;i<arg.length;i++){
            System.out.println(Arrays.toString(arg[i]));
        }
    }   
   
   
   
   
   
    public static boolean checking(ArrayList <int []> input, int [] arg) {
        boolean b=false;
        for(int i=0;i<input.size();i++){
            if(input.get(i)[0]==arg[0]&&input.get(i)[1]==arg[1]){
                b=true;
                break;
            }   
        }
          return b;
    }   
   
   

    public static String ToBinary(int input, int loci) {
        String b= Integer.toBinaryString(input);
        int holder=b.length();
       
        while(holder<loci){
            String h="0";
            b=h.concat(b);
            holder++;
        }
          return b;
    }   
   
    public static double [] Fitness(double  selection, int loci,int dimen) {
        double [] fitness=new double [dimen];
        for(int j=0;j<dimen;j++){
            String temp=ToBinary(j,loci);
            int holder=0;
            for(int k=0;k<loci;k++){
                if(temp.charAt(k)=='0'){holder++;}               
            }
            fitness[j]=selection*holder+0;
        }
        return fitness;
    }
   
   
    public static PrintWriter MakeFile(String input) {
        File file = new File(input);
        file.delete();
        File file1 = new File(input);
        FileWriter fw=null;
        try {
            fw = new FileWriter(file1, true);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            System.out.println(e);}
        BufferedWriter bw = new BufferedWriter(fw);
        PrintWriter out = new PrintWriter(bw, true);
        return out;
    }
   
    }

        // TODO Auto-generated method stub