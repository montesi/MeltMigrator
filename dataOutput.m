function Result=dataOutput(SpreadingRate_Half,ExtractionWidth,ExtractionDepth,SmoothingWidth,SeafloorCrust)

Result(1,1)=1;
Result(1,2)=SpreadingRate_Half/(3.17e-10);
Result(1,3)=ExtractionWidth;
Result(1,4)=ExtractionDepth;
Result(1,5)=SmoothingWidth;
Result(1,6)=SeafloorCrust.AverageWholdDomain;
Result(1,7)=SeafloorCrust.AverageRidgeDomain;
Result(1,8)=SeafloorCrust.AverageTransformDomain;
Result(1,9)=SeafloorCrust.AverageObliqueDomain;
    
return
    
    