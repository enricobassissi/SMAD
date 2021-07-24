function Population = int_pop_2RL_moo_AC_flex(GenomeLength, ~, options)
% INT_POP Function that creates an initial population satisfying bounds and
% integer constraints

totalPopulation = sum(options.PopulationSize);

%IntCon constraints
IntCon = [6,7,8,9,10];

range = options.PopInitRange;
lower = range(1,:);
span =  range(2,:) - lower;


Population = repmat(lower,totalPopulation,1 )+  ...
    repmat(span,totalPopulation,1) .* rand(totalPopulation, GenomeLength);

x = rand;
if x>=0.5
    Population(:,IntCon) = floor(Population(:, IntCon));
else
    Population(:,IntCon) = ceil(Population(:, IntCon));
end
Population = checkboundsIntGA(Population, range);