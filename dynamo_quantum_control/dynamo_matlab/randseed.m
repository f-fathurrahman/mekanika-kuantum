function randseed (seedvalue)

if ~exist('seedvalue','var')
    seedvalue = sum(100*clock);
end

%RandStream.setDefaultStream(RandStream('mt19937ar','seed', seedvalue));
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', seedvalue));