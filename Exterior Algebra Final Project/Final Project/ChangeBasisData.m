classdef (ConstructOnLoad) ChangeBasisData < event.EventData
   properties
       oldBasis
       newBasis
   end
   methods
       function eventData = ChangeBasisData(oldBasis,newBasis)
           eventData.oldBasis = oldBasis;
           eventData.newBasis = newBasis;
      end
   end
end