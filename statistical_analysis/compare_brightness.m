movieFiles={
    '20170724_094625', ...
    '20170724_100129', ...
    '20170724_101503', ...
    '20170724_102838', ...
    '20170724_104348', ...
    '20170725_092500', ...
    '20170725_093801', ...
    '20170725_095032', ...
    '20170725_100338', ...
    '20170725_101701', ...
    '20170726_103537', ...
    '20170726_104856', ...
    '20170726_110209', ...
    '20170726_111616', ...
    '20170726_113000', ...
    '20170727_094127', ...
    '20170727_095434', ...
    '20170727_100818', ...
    '20170727_102228', ...
    '20170727_103727', ...
    '20170728_095108'
    };

path = '~/cb/Data/';
time_cut = 30; % (cut first 30 secs)
fs = 5; % Hz
maxIter=10000;

numFiles = length(movieFiles);
%%

varlog=zeros(1,numFiles);

for m=1:numFiles
    filename = fullfile(path,[movieFiles{m} '.mat']);
    load(filename);
    
    fs=5; time_cut=30;
    yy = single(Y(141:1080-140,171:1440-170,1+time_cut*fs:end));
    clear Y;
    yy = reshape(yy,[],size(yy,3));
    yy=yy-min(yy);
    yy=yy./max(yy);
    
    varlog(m)=var(mean(yy));
    fprintf('Movie %d variance: %f\n',m,varlog(m));
end