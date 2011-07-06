#
#  numeric_utils.rb
#  ruby-genomics
#
#  Helper methods for determining numeric objects
#
#  Created by Timothy Palpant on 6/23/11.
#  Copyright 2011 UNC. All rights reserved.
#

class Object
  def integer?
    true if Integer(self) rescue false
  end
  
  def numeric?
    true if Float(self) rescue false
  end
end
