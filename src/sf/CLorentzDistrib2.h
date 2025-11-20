#ifndef C_LORENTZ_DISTRIB_H
#define C_LORENTZ_DISTRIB_H
#include <cmath>
#include "generatormt.h"

///Lorentz distribution of x:
///coefficient*halfWidth/( (x-centerValue)**2 + (halfWidth)**2)
///when the 1st constructor is used, the distribution is normalized on the [i_xMin; i_xMax] interval
///when the 2nd constructor is used, the distribution is normalized on the [-infinity; infinity]

class CLorentzDistrib
{
public:
	CLorentzDistrib(const double i_coeff,
        const double i_centerValue,
		const double i_halfWidth,
		const double i_xMin,
		const double i_xMax)
		:m_coeff( i_coeff ),
		 m_centerValue( i_centerValue ),
		 m_halfWidth( i_halfWidth ),
		 m_argMin( atan((i_xMin - i_centerValue)/m_halfWidth) ),
		 m_argDelta( atan((i_xMax - i_centerValue)/m_halfWidth) - m_argMin )
	{
	};

	CLorentzDistrib(const double i_coeff,
        const double i_centerValue,
		const double i_halfWidth)
		:m_coeff( i_coeff ),
		 m_centerValue( i_centerValue ),
		 m_halfWidth( i_halfWidth ),
		 m_argMin( -2.0*atan(1.0) ),
		 m_argDelta( 4.0*atan(1.0) )
	{
	};


	inline double eval(const double i_x) const
	{
		double arg( (m_centerValue - i_x)/m_halfWidth );

		return m_coeff/( ( 1.0 + arg*arg )*m_halfWidth*m_argDelta );
	};


	inline double generate() const
	{
		return m_centerValue + m_halfWidth*tan( m_argDelta*frandom() + m_argMin );
	};


private:
	const double m_coeff;
	const double m_centerValue;
	const double m_halfWidth;
	const double m_argMin;
	const double m_argDelta;
};

#endif
