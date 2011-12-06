#ifndef C_LORENTZ_DISTRIB_H
#define C_LORENTZ_DISTRIB_H

#include <cmath>

//Lorentz distribution of variable:
//halfWidth/( (x-centerValue)**2 + (halfWidth)**2)

class CLorentzDistrib
{
public:
	CLorentzDistrib(
		const double i_centerValue, 
		const double i_halfWidth)
		:
		 m_centerValue(i_centerValue),
		 m_halfWidth(i_halfWidth)
	{
	};

	inline double eval(const double i_x) const
	{
		double arg( (m_centerValue - i_x)/m_halfWidth );
		return 1.0/( ( 1.0 + arg*arg )*m_halfWidth );	
	};	


	inline double generate() const
	{	
		return m_centerValue + m_halfWidth*tan( M_PI*(frandom() - 0.5) );	
	};

	inline double generatePositiv() const
	{	
		return m_centerValue + m_halfWidth*tan( 0.5*M_PI*frandom() );	
	};

private:
	const double m_centerValue;
	const double m_halfWidth;
};


#endif
