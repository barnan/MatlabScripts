function retVal = ExtrPos3(left, middle, right, index)

    r2 = left + right - 2* middle;
    
    if abs(r2) > 1.0E-4     
        retVal = max(min(0.5*(left-right) /r2, 1.5), -1.5) + index;
    else
        retVal = index;        
    end
end


%from K.J. C# code:
%public static double ExtrPos3(double left, double middle, double right, int index)
%        {
%            double r2 = left + right - 2 * middle;
%            return Math.Abs(r2) > 1.0E-4D ? Math.Max(Math.Min(0.5 * (left - right) / r2, 1.50D), -1.50D) + index : index;
%        }