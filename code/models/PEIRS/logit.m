function out = logit(p)
    if  p < 0 | p > 1
      error('inverse logit only defined for numbers between 0 and 1');
    else
          out=log(p./(1-p));
    end
end

