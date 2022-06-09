
%% file paths
clear;
paths = {
'Z:\Alpha data\Neurobehavior\IC\ex110323\tank110323a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex110322\tank110322a\OurData-8'
'Z:\Alpha data\Neurobehavior\IC\ex110324\tank110324a\OurData-4'
'Z:\Alpha data\Neurobehavior\IC\ex110325\tank110325a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex110328\tank110328a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110329\tank110329a\OurData-8'
'Z:\Alpha data\Neurobehavior\IC\ex110330\tank110330a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110331\tank110331a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex110401\tank110401a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110407\tank110407a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110408\tank110408a\OurData-7'
%'Z:\Alpha data\Neurobehavior\IC\ex110411\tank110411a\OurData-6' DONT USE
'Z:\Alpha data\Neurobehavior\IC\ex110412\tank110412a\OurData-9'
%'Z:\Alpha data\Neurobehavior\IC\ex110413\tank110413a\OurData-7'DONT USE
'Z:\Alpha data\Neurobehavior\IC\ex110518\tank110518a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex110601\tank110601a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex110606\tank110606a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex100809\tank100809a\OurData-3'
'Z:\Alpha data\Neurobehavior\IC\ex100917\tank100917a\OurData-2'
'Z:\Alpha data\Neurobehavior\IC\ex101021\tank101021a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex101022\tank101022a\OurData-6'
'Z:\Alpha data\Neurobehavior\IC\ex101025\tank101025a\OurData-7'
'Z:\Alpha data\Neurobehavior\IC\ex101026\tank101026a\OurData-4'
'Z:\Alpha data\Neurobehavior\IC\ex101027\tank101027a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex101028\tank101028a\OurData-7'
'Z:\Alpha data\Neurobehavior\IC\ex101029\tank101029a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex101102\tank101102a\OurData-4'
'Z:\Alpha data\Neurobehavior\IC\ex101103\tank101103a\OurData-9'
'Z:\Alpha data\Neurobehavior\IC\ex101105\tank101105a\OurData-7'
'Z:\Alpha data\Neurobehavior\IC\ex101117\tank101117a\OurData-10'
'Z:\Alpha data\Neurobehavior\IC\ex110120\tank110120a\OurData-10'
'Z:\Alpha data\Neurobehavior\IC\ex110121\tank110121a\OurData-10'
'Z:\Alpha data\Neurobehavior\IC\ex110125\tank110125a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110126\tank110126a\OurData-5'
'Z:\Alpha data\Neurobehavior\IC\ex110127\tank110127a\OurData-4'
'Z:\Alpha data\Neurobehavior\IC\ex110128\tank110128a\OurData-4'
'Z:\Alpha data\Neurobehavior\IC\ex110209\tank110209a\OurData-4'
%'Z:\Alpha data\Neurobehavior\CN\ex110808\tank110808a\OurData-6'DONT USE
'Z:\Alpha data\Neurobehavior\CN\ex110811\tank110811a\OurData-4'
%'Z:\Alpha data\Neurobehavior\CN\ex110815\tank110815a\OurData-5'DONT USE
'Z:\Alpha data\Neurobehavior\CN\ex110826\tank110826a\OurData-8'
'Z:\Alpha data\Neurobehavior\CN\ex110902\tank110902a\OurData-5'
'Z:\Alpha data\Neurobehavior\CN\ex110907\tank110907a\OurData-5'
'Z:\Alpha data\Neurobehavior\CN\ex110909\tank110909a\OurData-6'
'Z:\Alpha data\Neurobehavior\CN\ex110919\tank110919a\OurData-5'
'Z:\Alpha data\Neurobehavior\CN\ex110922\tank110922a\OurData-8'
'Z:\Alpha data\Neurobehavior\CN\ex110923\tank110923a\OurData-10'
'Z:\Alpha data\Neurobehavior\CN\ex111019\tank111019a\OurData-8'
'Z:\Alpha data\Neurobehavior\CN\ex111021\tank111021a\OurData-10'
'Z:\Alpha data\Neurobehavior\CN\ex111024\tank111024a\OurData-7'

% 'Z:\Bravo data\Neurobehavior\IC\ex100304\tank100304b\OurData-1'
% 'Z:\Bravo data\Neurobehavior\IC\ex100305\tank100305b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\IC\ex100308\tank100308b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\IC\ex100310\tank100310b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\IC\ex100310\tank100310b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100325\tank100325b\OurData-10'
% 'Z:\Bravo data\Neurobehavior\IC\ex100325\tank100325b\OurData-9'
% 'Z:\Bravo data\Neurobehavior\IC\ex100325\tank100325b\OurData-3'
% 'Z:\Bravo data\Neurobehavior\IC\ex100326\tank100326b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100330\tank100330b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\IC\ex100707\tank100707b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\IC\ex100708\tank100708b\OurData-6'
% 'Z:\Bravo data\Neurobehavior\IC\ex100712\tank100712b\OurData-6'
% 'Z:\Bravo data\Neurobehavior\IC\ex100713\tank100713b\OurData-10'
% 'Z:\Bravo data\Neurobehavior\IC\ex100714\tank100714b\OurData-14'
% 'Z:\Bravo data\Neurobehavior\IC\ex100720\tank100720b\OurData-9'
% 'Z:\Bravo data\Neurobehavior\IC\ex100720\tank100720b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\IC\ex100721\tank100721b\OurData-3'
% 'Z:\Bravo data\Neurobehavior\IC\ex100723\tank100723b\OurData-1'
% 'Z:\Bravo data\Neurobehavior\IC\ex100804\tank100804b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100804\tank100804b\OurData-6'
% 'Z:\Bravo data\Neurobehavior\IC\ex100812\tank100812b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100812\tank100812b\OurData-7'
% 'Z:\Bravo data\Neurobehavior\IC\ex100813\tank100813b\OurData-3'
% 'Z:\Bravo data\Neurobehavior\IC\ex100818\tank100818b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100915\tank100915b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100916\tank100916b\OurData-2'
% 'Z:\Bravo data\Neurobehavior\IC\ex100917\tank100917b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\IC\ex101015\tank101015b\OurData-8'
% 'Z:\Bravo data\Neurobehavior\IC\ex101020\tank101020b\OurData-7'
% 'Z:\Bravo data\Neurobehavior\IC\ex101020\tank101020b\OurData-26'
% 'Z:\Bravo data\Neurobehavior\IC\ex101022\tank101022b\OurData-8'
% 'Z:\Bravo data\Neurobehavior\IC\ex101027\tank101027b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\IC\ex101027\tank101027b\OurData-14'
% 'Z:\Bravo data\Neurobehavior\IC\ex101101\tank101101b\OurData-7'
% 'Z:\Bravo data\Neurobehavior\CN\ex111003\tank111003b\OurData-7'
% %'Z:\Bravo data\Neurobehavior\CN\ex111006\tank111006b\OurData-6'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex111007\tank111007b\OurData-5'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex111010\tank111010b\OurData-7'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex111011\tank111011b\OurData-6'
% %'Z:\Bravo data\Neurobehavior\CN\ex110602\tank110602b\OurData-6'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex110421\tank110421b\OurData-7'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex110420\tank110420b\OurData-8'
% 'Z:\Bravo data\Neurobehavior\CN\ex110418\tank110418b\OurData-13'
% 'Z:\Bravo data\Neurobehavior\CN\ex110413\tank110413b\OurData-9'
% 'Z:\Bravo data\Neurobehavior\CN\ex110412\tank110412b\OurData-12'
% 'Z:\Bravo data\Neurobehavior\CN\ex110411\tank110411b\OurData-13'
% 'Z:\Bravo data\Neurobehavior\CN\ex110408\tank110408b\OurData-6'
% 'Z:\Bravo data\Neurobehavior\CN\ex110407\tank110407b\OurData-6'
% %'Z:\Bravo data\Neurobehavior\CN\ex110406\tank110406b\OurData-7'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex110405\tank110405b\OurData-11'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex110404\tank110404b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\CN\ex110404\tank110404b\OurData-7'
% %'Z:\Bravo data\Neurobehavior\CN\ex110401\tank110401b\OurData-12'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex110331\tank110331b\OurData-8'
% % 'Z:\Bravo data\Neurobehavior\CN\ex110606\tank110606b\OurData-16'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex110608\tank110608b\OurData-7'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex110610\tank110610b\OurData-6'
% %'Z:\Bravo data\Neurobehavior\CN\ex110614\tank110614b\OurData-6'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex110729\tank110729b\OurData-8'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex110805\tank110805b\OurData-8'
% 'Z:\Bravo data\Neurobehavior\CN\ex110808\tank110808b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\CN\ex110829\tank110829b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\CN\ex110830\tank110830b\OurData-8'
% %'Z:\Bravo data\Neurobehavior\CN\ex110906\tank110906b\OurData-5'DONT USE
% %'Z:\Bravo data\Neurobehavior\CN\ex111012\tank111012b\OurData-5'DONT USE
% 'Z:\Bravo data\Neurobehavior\CN\ex111102\tank111102b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\CN\ex110330\tank110330b\OurData-4'
% 'Z:\Bravo data\Neurobehavior\CN\ex110329\tank110329b\OurData-9'
% 'Z:\Bravo data\Neurobehavior\CN\ex110328\tank110328b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\CN\ex110323\tank110323b\OurData-5'
% 'Z:\Bravo data\Neurobehavior\CN\ex101118\tank101118b\OurData-5'

% 'Z:\Charlie data\Neurobehavior\IC\ex121016\tank121016c\OurData-4'
% 'Z:\Charlie data\Neurobehavior\IC\ex121109\tank121109c\OurData-9'
% 'Z:\Charlie data\Neurobehavior\IC\ex121206\tank121206c\OurData-5'
% 'Z:\Charlie data\Neurobehavior\IC\ex121207\tank121207c\OurData-8'
% 'Z:\Charlie data\Neurobehavior\IC\ex130114\tank130114c\OurData-6'
% %'Z:\Charlie data\Neurobehavior\IC\ex130118\tank130118c\OurData-5'DONTUSE
% 'Z:\Charlie data\Neurobehavior\IC\ex130122\tank130122c\OurData-7'
% 'Z:\Charlie data\Neurobehavior\IC\ex130201\tank130201c\OurData-8'
% 'Z:\Charlie data\Neurobehavior\IC\ex130214\tank130214c\OurData-6'
% 'Z:\Charlie data\Neurobehavior\IC\ex130215\tank130215c\OurData-9'
% 'Z:\Charlie data\Neurobehavior\IC\ex130308\tank130308c\OurData-13'
% 'Z:\Charlie data\Neurobehavior\IC\ex130417\tank130417c\OurData-13'
% 'Z:\Charlie data\Neurobehavior\IC\ex130418\tank130418c\OurData-7'
% 'Z:\Charlie data\Neurobehavior\IC\ex130625\tank130625c\OurData-7'
% 'Z:\Charlie data\Neurobehavior\IC\ex130701\tank130701c\OurData-5'
% 'Z:\Charlie data\Neurobehavior\IC\ex130828\tank130828c\OurData-15'
% 'Z:\Charlie data\Neurobehavior\IC\ex130923\tank130923c\OurData-8'
% 'Z:\Charlie data\Neurobehavior\IC\ex140314\tank140314c\OurData-10'
% 'Z:\Charlie data\Neurobehavior\IC\ex140917\tank140917c\OurData-12'
% 'Z:\Charlie data\Neurobehavior\CN\ex141111\tank141111c\OurData-4'

% 'Z:\Delta data\Neurobehavior\IC\ex120822\tank120822d\OurData-10'
% 'Z:\Delta data\Neurobehavior\IC\ex120823\tank120823d\OurData-10'
% 'Z:\Delta data\Neurobehavior\IC\ex120827\tank120827d\OurData-6'
% 'Z:\Delta data\Neurobehavior\IC\ex120907\tank120907d\OurData-7'
% 'Z:\Delta data\Neurobehavior\IC\ex120911\tank120911d\OurData-6'
% 'Z:\Delta data\Neurobehavior\IC\ex120927\tank120927d\OurData-8'
% 'Z:\Delta data\Neurobehavior\IC\ex121004\tank121004d\OurData-9'
% 'Z:\Delta data\Neurobehavior\IC\ex121012\tank121012d\OurData-7'
% 'Z:\Delta data\Neurobehavior\IC\ex121016\tank121016d\OurData-7'
% 'Z:\Delta data\Neurobehavior\IC\ex121017\tank121017d\OurData-9'
% 'Z:\Delta data\Neurobehavior\IC\ex121030\tank121030d\OurData-8'
% 'Z:\Delta data\Neurobehavior\IC\ex121101\tank121101d\OurData-6'
% 'Z:\Delta data\Neurobehavior\IC\ex121108\tank121108d\OurData-11'
% 'Z:\Delta data\Neurobehavior\IC\ex121115\tank121115d\OurData-10'
% 'Z:\Delta data\Neurobehavior\IC\ex121116\tank121116d\OurData-6'
% 'Z:\Delta data\Neurobehavior\IC\ex121212\tank121212d\OurData-9'
% 'Z:\Delta data\Neurobehavior\IC\ex130114\tank130114d\OurData-7'
% 'Z:\Delta data\Neurobehavior\IC\ex130116\tank130116d\OurData-7'
% 'Z:\Delta data\Neurobehavior\IC\ex121113\tank121113d\OurData-11'
% 'Z:\Delta data\Neurobehavior\IC\ex121206\tank121206d\OurData-10'
% 'Z:\Delta data\Neurobehavior\CN\ex140414\tank140414d\OurData-11'
% 'Z:\Delta data\Neurobehavior\CN\ex140505\tank140505d\OurData-14'
% 'Z:\Delta data\Neurobehavior\CN\ex140409\tank140409d\OurData-11'
% 'Z:\Delta data\Neurobehavior\CN\ex140429\tank140429d\OurData-12'
% 'Z:\Delta data\Neurobehavior\CN\ex140317\tank140317d\OurData-12'
% 'Z:\Delta data\Neurobehavior\CN\ex140318\tank140318d\OurData-11'
% 'Z:\Delta data\Neurobehavior\CN\ex140319\tank140319d\OurData-17'
% 'Z:\Delta data\Neurobehavior\CN\ex140325\tank140325d\OurData-9'
% 'Z:\Delta data\Neurobehavior\CN\ex140723\tank140723d\OurData-17'
% 'Z:\Delta data\Neurobehavior\CN\ex140929\tank140929d\OurData-9'
% 'Z:\Delta data\Neurobehavior\CN\ex140807\tank140807d\OurData-8'


        };
%% import all of the above data into single cell variable 

%Num_monkeys = 3; %Figure out how to make a 3-d cell array, third dimension is Number of monkeys (3)
% Each row is new tank/block combo  
neuron_data = cell(length(paths),1);

parfor imp = 1:length(paths)
    neuron_data{imp,1} = TDTbin2mat(paths{imp,1});
end

%% Loop trough each of the duration values and Calc Grand_CP for all durations 


dur = 12.5:12.5:200;

CP = cell(length(neuron_data),length(dur));
Grand_CP = cell(length(neuron_data),length(dur)); 
p_values = cell(length(neuron_data),length(dur));

for n_ct = 1:length(neuron_data)
    for dur_ct = 1:length(dur)
           [CP{n_ct,dur_ct},Grand_CP{n_ct,dur_ct},p_values{n_ct,dur_ct}] = CP_Calc(neuron_data{n_ct,1},dur(dur_ct));
    end
end

    
    
    
    
    
    