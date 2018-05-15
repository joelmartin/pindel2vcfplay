#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;


/* 'Parameter' stores an individual parameter; it is set by the command line parameters, and used by the program. */
class Parameter {
public:
	bool isRequired() const {
		return d_required;
	}
	void describe() const;
	string getDescription() const {
		return d_description;
	}
	string getShortName() const {
		return d_shortName;
	}
	string getLongName() const {
		return d_longName;
	}
	bool hasName( const string& name ) const {
		return ( d_shortName.compare(name) == 0 || d_longName.compare( name ) == 0 );
	}
	bool isSet() const {
		return d_isSet;
	}
	virtual void setValue(const string& value ) {
		cout << "WHAT!" << endl ;
	};
	virtual void setValue(const int value) {};
	virtual void setValue(const double value) {};
	virtual void setValue(const bool value) {};
	virtual int getIValue() const {};
	virtual bool getBValue() const {};
	virtual string getSValue() const {};
	virtual double getFValue() const {};

	virtual bool isUnary() const {
		return false;
	}

	Parameter( const string& shortName, const string& longName, const string& description, const bool required);
protected:
	void set() {
		d_isSet = true;
		//cout << "setting " << d_shortName << endl;
	}


private:
	bool d_required;
	string d_shortName;
	string d_longName;
	string d_description;
	bool d_isSet;
};

class IntParameter: public Parameter {
public:
	IntParameter( int* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const int value );

	virtual int getIValue() const {
		return *d_data_ptr;
	}
	virtual void setValue(const string& value);
	virtual void setValue(const int value);
private:
	int* d_data_ptr;
};

class BoolParameter: public Parameter {
public:
	BoolParameter( bool* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const bool value );

	virtual bool getBValue() const {
		return *d_data_ptr;
	}
	virtual void setValue(const string& value);
	virtual void setValue(const bool value);
	virtual bool isUnary() const {
		return true;
	}
private:
	bool* d_data_ptr;
};

class FloatParameter: public Parameter {
public:
	FloatParameter( double* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const double value );

	double getFValue() const {
		return *d_data_ptr;
	}
	void setValue(const string& value);
	void setValue(const double value);
private:
	double* d_data_ptr;
};

class StringParameter: public Parameter {
public:
	StringParameter( string* par_ptr, const string& shortName, const string& longName, const string& description, const bool required, const string& value );

	string getSValue() const {
		return *d_data_ptr;
	}
	void setValue(const string& value);
private:
	string* d_data_ptr;
};

#endif
