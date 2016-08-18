#include "st_ed.h"

st_ed::st_ed(int pos, char ids)
{
	_position = pos;
	_ids = ids;
}

st_ed::st_ed()
{
	_position = 0;
	_ids = 's';
}

st_ed::~st_ed()
{
}

bool st_ed::PointerCompare::operator()(const st_ed* l, const st_ed* r) 
{
	return *l < *r;
}

bool st_ed::operator<(const st_ed &other) const
{
	if(_position < other._position)
	{
		return 1;
	}
	else
		return 0;
}
bool st_ed::operator>(const st_ed &other) const
{
	if(_position > other._position)
	{
		return 1;
	}
	else
		return 0;
}
bool st_ed::operator==(const st_ed &other) const
{
	if(_position == other._position)
	{
		return 1;
	}
	else
		return 0;
}
bool st_ed::operator!=(const st_ed &other) const
{
	if(_position != other._position)
	{
		return 1;
	}
	else
		return 0;
}