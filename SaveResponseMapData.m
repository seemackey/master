function []=SaveResponseMapData(x,y,baseline,Monkey,Loc,TrialNum,DataTank)
   TrialNum=num2str(TrialNum);
   Dirrectory= ['C:\Users\Kimaya\Desktop\ResponseMapData\',Monkey,'_',Loc,'_',DataTank,'_TrialNum',TrialNum,'_ResponseMap.txt'];
   lengthResp=length(x); 
   a=fopen(Dirrectory, 'w');
   fprintf(a, '%d', lengthResp);
   fprintf(a, '\r\n');
   fprintf(a, '%d\t', x);
   fprintf(a, '\r\n');
   fprintf(a, '%d\t', y);
   fprintf(a, '\r\n');
   fprintf(a, '%d\t', baseline);
   fclose(a);
end   
   
   
