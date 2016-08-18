#ifndef ST_ED_H
#define ST_ED_H

class st_ed
{
public:
	st_ed(int pos, char ids);
	st_ed();
	~st_ed();
	int _position;
	int _dest;
	char _ids;
	
	bool operator<(const st_ed &other) const;
	bool operator>(const st_ed &other) const;
	bool operator==(const st_ed &other) const;
	bool operator!=(const st_ed &other) const;
	
	struct PointerCompare {
      bool operator()(const st_ed* l, const st_ed* r);
	};
};

#endif // ST_ED_H
