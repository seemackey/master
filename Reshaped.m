function [new_tlvlsC,new_tlvlsI] = Reshaped(tlvlsCorrect,tlvlsIncorrect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Creates Tone level Data with >=3 trials at the same tone level(s) %%%%%%%%%%%%%%%%%%%%%
%%%%%%   Outputs the Reshaped Version of tlvlsCorrect & tlvlsIncorrect
%%%%%%     

%%%%%%%%%%%%%%%%%%%%%%%%%%Jackson Mayfield 08/30/21%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[C_rows,C_columns] = size(tlvlsCorrect); %Initial Size of data matrices for indexing
[I_rows,I_columns] = size(tlvlsIncorrect);  

%%%%%%%%%%%%%% Correct Tone level Data Reshaping 
counter_col = zeros; 
counter = 0;
for i = 1:C_columns    %Counts all the numeric values only in each of the data column
    for j = 1:C_rows
        if (isnan(tlvlsCorrect(j,i)) == 0) && (tlvlsCorrect(j,i) ~= 0)
            counter = counter + 1;   
        end
         counter_col(i) = counter;
    end
    counter = 0; 
end

i = 1; 
j = 1;
del_col = zeros;
while i <= length(counter_col)   %Removes all tone level data sets with less than 3 trials
    if counter_col(i) < 4
            del_col(j) = i; 
            j = j + 1; 
    end
    i = i + 1; 
end

%%%%%%%%%%%%%%%% Incorrect Tone Level Data Reshaping
counter_col2 = zeros; 
counter2 = 0; 
for i = 1:I_columns
    for j = 1:I_rows
        if (isnan(tlvlsIncorrect(j,i)) == 0) && (tlvlsIncorrect(j,i) ~= 0)
            counter2 = counter2 + 1;  
        end
        counter_col2(i) = counter2;
    end
    counter2 =0;
end

i = 1; 
j = 1;
del_col2 = zeros;
while i <= length(counter_col2)  
    if counter_col2(i) < 4
            del_col2(j) = i; 
            j = j + 1; 
    end
    i = i + 1; 
end

new_tlvlsC = tlvlsCorrect; %Just so we can rename for future use in the code
new_tlvlsI = tlvlsIncorrect;

if del_col > 0
new_tlvlsC(:,del_col) = []; %New tone levels correct with >= 3 trials 
end

if del_col2 > 0
new_tlvlsI(:,del_col2) = []; %New tone levels incorrect with >= 3 trials 
end

%If the first row (tone level) is the same for Correct and Incorrect, then add to histogram/ 
%if not you cant add both C and I to same tone level Histogram 
k = 1;
hist_trial_Col = zeros;
for i = 1:length(new_tlvlsC(1,:))
    for j = 1:length(new_tlvlsI(1,:))
        if  new_tlvlsC(1,i) == new_tlvlsI(1,j)
            hist_trial_Col(k) = new_tlvlsC(1,i); %gets the value(s) for the tone levels that are the same on Both C and I
            k = k + 1;
        end  
    end
   
end
if hist_trial_Col == 0
    
    disp('Correct and Incorrect Responses do not fall under the same tone levels. Remove tank from Analysis.')
    return
end

        q = 1;
        columnC = zeros;
        for i = 1:length(new_tlvlsC(1,:))
            for j = 1:length(hist_trial_Col(1,:))
                if new_tlvlsC(1,i) == hist_trial_Col(1,j)
                    columnC(q) = i;
                    q = q + 1;
                end
            end
        end
        q = 1;
        columnI = zeros; 
        for i = 1:length(new_tlvlsI(1,:))
            for j = 1:length(hist_trial_Col(1,:))
                if new_tlvlsI(1,i) == hist_trial_Col(1,j)
                    columnI(q) = i;
                    q = q + 1;
                end
            end
        end


new_tlvlsC = new_tlvlsC(:,columnC);
new_tlvlsI = new_tlvlsI(:,columnI);

%% new edits - CM
% NaN indices 
[Cnanrow,Cnancol]=find(isnan(new_tlvlsC));
[Inanrow,Cnancol]=find(isnan(new_tlvlsI));

% delete anything after the NaNs b/c this fxn inserts zeros sometimes
if length(new_tlvlsI)>length(new_tlvlsC)
   new_tlvlsI(max(Cnanrow)+1:end,:)=[];
end

if length(new_tlvlsC)>length(new_tlvlsI)
    new_tlvlsC(max(Inanrow)+1:end,:)=[];
end


%new_tlvlsC = rmmissing(new_tlvlsC) %removes NaN values from single column or trial
%new_tlvlsI = rmmissing(new_tlvlsI) 


end
