%{
#include "BlowBotl.h"
%}

namespace stk {
    class BlowBotl : public Instrmnt
    {
    public:

        BlowBotl( void );
        ~BlowBotl( void );

        void clear( void );
        void setFrequency( StkFloat frequency );
        void startBlowing( StkFloat amplitude, StkFloat rate );
        void stopBlowing( StkFloat rate );
        void noteOn( StkFloat frequency, StkFloat amplitude );
        void noteOff( StkFloat amplitude );
        void controlChange( int number, StkFloat value );
        StkFloat tick( unsigned int channel = 0 );

        StkFrames& tick( StkFrames& frames, unsigned int channel = 0 );
    };
}
