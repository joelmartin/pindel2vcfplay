#include <string>
#include <iostream>
#include "p2v_parameters.h"

IntParameter::IntParameter( int* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const int value ) :
   Parameter( shortName, longName, description, required ), d_data_ptr( par_ptr )
{
   *d_data_ptr = value;
}

void IntParameter::setValue(const string& value)
{
   setValue( atoi( value.c_str() ));
}

void IntParameter::setValue(const int value)
{
   *d_data_ptr = value;
   set();
}

BoolParameter::BoolParameter( bool* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const bool value ) :
   Parameter( shortName, longName, description, required )
{
   d_data_ptr = par_ptr;
   *d_data_ptr = value;
}

void BoolParameter::setValue(const string& value)
{
   char firstChar = tolower( value[0] );
   setValue( ( firstChar == 'f' || firstChar == '0' ) ? false : true );
}

void BoolParameter::setValue(const bool value)
{
   *d_data_ptr = value;
   set();
}

FloatParameter::FloatParameter( double* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const double value ) :
   Parameter( shortName, longName, description, required )
{
   d_data_ptr = par_ptr;
   *d_data_ptr = value;
}

void FloatParameter::setValue(const string& value)
{
   setValue( atof( value.c_str() ));
}

void FloatParameter::setValue(const double value)
{
   *d_data_ptr = value;
   set();
}

StringParameter::StringParameter( string* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const string& value ) :
   Parameter( shortName, longName, description, required )
{
   d_data_ptr = par_ptr;
   *d_data_ptr = value;
}

void StringParameter::setValue(const string& value)
{
   *d_data_ptr = value;
   set();
}

Parameter::Parameter( const string& shortName, const string& longName, const string& description, const bool required )
{
   d_required = required;
   d_shortName = shortName;
   d_longName = longName;
   d_description = description;
   d_isSet = false;
}

void Parameter::describe() const
{
   cout << d_shortName << "/" << d_longName << "  " << d_description;
   if ( d_required ) {
      cout << ": required parameter" ;
   }
   cout << endl;
}

