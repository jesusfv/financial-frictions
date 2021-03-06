% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


% start with 1, then override to 0 if a check fails
use_parallel = 1;

% first check if the parallel toolbox is installed
try 
   ver('distcomp')
catch
   disp('Parallel Computing Toolbox not installed')
   use_parallel = 0;
end
disp(' ');


% if there's less than 12 GB of RAM, better not to use the parallel toolbox even if it's available
% the "memory" function is only available on Windows, so if you're on Linux or Mac and have <12GB of RAM please use the manual override at the end

if ispc
    [userview,systemview] = memory;
    if systemview.PhysicalMemory.Total < 12*1024*1024*1024
        disp('Total memory <12GB');
        use_parallel = 0;
    end
end


% manual override
% use_parallel = 0;



