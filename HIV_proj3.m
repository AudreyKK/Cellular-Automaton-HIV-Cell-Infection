% =========================================================================
%  Final Project - Modeling HIV Via CA
%  PART III
%  
%  THE ONLY DIFFERENCE in this code is that probRespond is being used as
%  probRespond_func, a negatively linear function.  I thought it would be
%  easier to have a separate instance of the code rather than having
%  different commented-out sections of project part II with complicated
%  instructions on what to do.
%
%  Audrey Keefe
% =========================================================================

% =====
% 1a
% =====
DEAD = 0;
HEALTHY = 1;
INFECTED_A1 = 2;
INFECTED_A2 = 3;
EMPTY = -1; % for the external grid

% Initial probability of infected_A1 cells
probHIV = 0.05;

% probability that a dead cell will be replaced by a healthy cell at the
% next time step
probReplace = 0.99;

% probability that a new healthy cell will be replaced by an infected_A1
% cell
probInfect = 10e-5;

responseTime = 4; % weeks until immune response

% Grid dimensions
m = 70;
n = 70;

% Creates a grid of 1s noting location of INFECTED_A1 cells
init_infected = (rand(m,n) < probHIV);

% Creates a grid of 1s noting location of HEALTHY cells
init_healthy = ones(m,n) - init_infected;

% Creates the initial state grid
cell_grid = init_healthy * HEALTHY + init_infected * INFECTED_A1;

% Create the external grid, with absorbing boundaries

extGrid = EMPTY * ones(m + 2, n + 2);

extGrid(2:m+1,2:n+1) = cell_grid;

% A number 0 <= rankLevel <= 8 (the number of neighbors) that indicates
% the effectiveness of a drug therapy (project 2)
rankLevel = 1;

% The probability that a patient will respond to a drug (project 2)
k = 1;
r = -0.01
probRespond_func(k) = r*k;

% =========================================================================
%  RULES
% =========================================================================
%
% 1) During drug therapy a HEALTHY cell with rankLevel number or more of
% INFECTED_A1 nighbors becomes infected_A1 with a probability of 
%   (1 - probRespond) * rankLevel/8), where probRespond is a
%   response-to-therapy related probability.
%
% 2) A 'healthy' cell with 'n' number of 'infected_A2' neighbors where
% 3 <= n <= 8, becomes 'infected_A1' because 'infected_A2' cells with
% concentration above a certain threshold can contaminate a 'healthy' cell.
%
% 3) All other 'healthy' cells remain healthy
%
% 4) An 'infected_A1' cell becomes an 'infected_A2' cell after
% 'responseTime', the number of time steps for the immune system to
% generate a response to kill an 'infected_A1' cell.
%
% 5) An 'infected_A2' cell becomes a 'dead' cell.
%
% 6) A 'dead' cell becomes a 'healthy' cell at the next time step with a
% probability of 'probReplace' because the immune system has a great
% ability to recover from an infection's immunosuppressant.
%
% 7) A new, 'healthy' cell may be replaced by an 'infected_A1' cell with a
% probability of 'probInfect' because new infected cells can come into the
% system.



cell_gridList{1} = cell_grid;
extGridList{1} = extGrid;

% =====
% 1b : plot numbers of healthy, infected, and dead cells versus time from 0
% through 12 weeks and then from 0 through 12 years.
% =====

% There are commented-out sections of the code.  These parts should be run
% when you are trying to test the system with a probRespond that declines
% over time.  In other words, when the probability of of a patient
% responding positively to a drug declines over time.
% This part starts only when the drug regime begins.  The probability of
% response declines every 26 weeks so that the person running the
% simulation can see the change in the system when using interval steps of
% 13 weeks.  I have arbitrarily chosen 13 weeks to be the ideal step
% because it allows the simulator to see the full scope of the changes
% without stepping page by page.

numIterations = 626; % week(s) 626 weeks for 12 years

for k=2:numIterations + 1
    cell_grid = cell_gridList{k-1};
    extGrid = extGridList{k-1};
    
    new_cellGrid = zeros(m,n);
    
    probRespond = probRespond_func(k-1);
    
    for j=2:n+1
        for i=2:m+1
            site = extGrid(i,j);
            N = extGrid(i - 1,j);
            NE = extGrid(i - 1,j + 1);
            NW = extGrid(i - 1,j - 1);
            E = extGrid(i,j + 1);
            S = extGrid(i + 1,j);
            SE = extGrid(i + 1, j + 1);
            SW = extGrid(i + 1, j - 1);
            W = extGrid(i,j - 1);
           
            A2_grid = [NW,N,NE;W,0,E;SW,S,SE];
            A2_count = sum(sum(A2_grid == INFECTED_A2));
            
            A1_grid = [NW,N,NE;W,0,E;SW,S,SE];
            A1_count = sum(sum(A1_grid == INFECTED_A1));
            
            if site == DEAD
                if rand < probReplace
                    new_cellGrid(i-1,j-1) = HEALTHY;
                else
                    new_cellGrid(i-1,j-1) = DEAD;
                end
            elseif site == HEALTHY
                if (N == INFECTED_A1 || NE == INFECTED_A1 || ...
                        NW == INFECTED_A1 || E == INFECTED_A1 ...
                        || W == INFECTED_A1 || SW == INFECTED_A1 ...
                        || S == INFECTED_A1 || SE ==INFECTED_A1)
                    % drug regime starts at week 300
                    if (k-1 >= 300) && (A1_count >= rankLevel )
                        % respondA
                        % Every 26 weeks, probRespond declines by 10%
                        if mod(k,26) == 0
                            if probRespond < 0
                                probRespond = 0;
                            end
                        end
                        %
                        if rand < (1 - probRespond) * rankLevel/8
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                        else
                        new_cellGrid(i-1,j-1) = HEALTHY;
                        end
                    else
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                    end   
                elseif A2_count >= 3
                    % respondA
                    % drug regime starts at week 300
                    if (k-1 >= 300) && (A1_count >= rankLevel)
                        % Every 26 weeks, probRespond declines by 10%
                        if mod(k,26) == 0
                            if probRespond < 0
                                probRespond = 0;
                            end
                        end
                        %
                        if rand < (1-probRespond) * rankLevel/8
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                        else
                            new_cellGrid(i-1,j-1) = HEALTHY;
                        end
                    else
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                    end
                elseif rand < probInfect
                    if (k-1 >= 300) && (A1_count >= rankLevel)
                        % respondA
                        % Every 26 weeks, probRespond declines by 10%
                        if mod(k,26) == 0
                            if probRespond < 0
                                probRespond = 0;
                            end
                        end
                        % 
                        if rand < (1-probRespond) * rankLevel/8
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                        else
                        new_cellGrid(i-1,j-1) = HEALTHY;
                        end
                    else
                        new_cellGrid(i-1,j-1) = INFECTED_A1;
                    end
                else
                    new_cellGrid(i-1,j-1) = HEALTHY;
                end
            elseif site == INFECTED_A1
                if mod(k,responseTime) == 0
                    new_cellGrid(i-1,j-1) = INFECTED_A2;
                else
                    new_cellGrid(i-1,j-1) = INFECTED_A1;
                end
            elseif site == INFECTED_A2
                new_cellGrid(i-1,j-1) = DEAD;
            end
            
        end
    end
    
    probRespond_func(k) = probRespond_func(k-1);
    cell_grid = new_cellGrid;
    extGrid(2:m+1,2:n+1) = cell_grid;
    cell_gridList{k} = cell_grid;
    extGridList{k} = extGrid;
end

avg_healthy = sum(sum(cell_gridList{1} == HEALTHY));
avg_A1 = sum(sum(cell_gridList{1} == INFECTED_A1));
avg_A2 = sum(sum(cell_gridList{1} == INFECTED_A2));

for i=2:length(extGridList)
    avg_healthy = avg_healthy + sum(sum(cell_gridList{i} == HEALTHY));
    avg_A1 = avg_A1 + sum(sum(cell_gridList{i} == INFECTED_A1));
    avg_A2 = avg_A2 + sum(sum(cell_gridList{i} == INFECTED_A2));

end

% average amount of healthy cells overall
final_avg_healthy = avg_healthy/length(extGridList);

% average amount of INFECTED_A1 cells overall
final_avg_A1 = avg_A1/length(extGridList);

% average amount of INFECTED_A2 cells overall
final_avg_A2 = avg_A2/length(extGridList);

set(groot,'DefaultFigureColormap',jet(64));
% Dark red    -> INFECTED_A2
% Orange      -> INFECTED_A1
% Light Green -> HEALTHY
% Light Blue  -> DEAD

upper = 3;
lower = 0;

interval = 13; % week(s). Interval is 1 for the 12 week period and 13 for
               % the 12 year period.

for i=1:interval:length(extGridList)
    
    data = extGridList{i};
    
     cmax = max(data(:));
    if cmax < upper
        cmax = upper;
    end
    cmin = min(data(:));
    if cmin > lower
        cmin = lower;
    end
    
    imagesc(data);
    caxis([cmin cmax]);
    colorbar;
    title(sprintf('Frame: %d',i));
    hold;
    axis equal; axis tight; axis xy;
    
    fprintf('Waiting for any key to be pressed.\n');
    w = waitforbuttonpress;
    
end