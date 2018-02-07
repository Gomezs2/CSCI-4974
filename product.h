#include <string>
#include <iostream>
#include <vector>

class Product {
public:
	Product(): id_(-1), asin_(""), title_(""), group_(""), avgRating_(-1) {}
	Product(int id): id_(id), asin_(""), title_(""), group_(""), avgRating_(-1) {}
	Product(int id, const std::string& asin, const std::string& title, const std::string& group, float avg_rating):
		id_(id), asin_(asin), title_(title), group_(group), avgRating_(avg_rating) {}

	int getId() const {return id_;}
	float getRating() const {return avgRating_;}
	const std::string& getTitle() const {return title_;}
	const std::string& getAsin() const {return asin_;}
	const std::string& getGroup() const {return group_;}
	const std::vector<std::string>& getSProducts() const {return Simproducts_;}

	void addSProduct(const std::string& item) {Simproducts_.push_back(item);}
	void addRating(float avg_rating) 	{avgRating_ = avg_rating;}
	void addTitle(const std::string& title)  {title_ = title;}
	void addAsin(const std::string& asin)   {asin_ = asin;}
	void addGroup(const std::string& group)  {group_ = group;}
	void print() const;

private:
	int id_;
	float avgRating_;
	std::string title_;
	std::string asin_;
	std::string group_;
	std::vector<std::string> Simproducts_;
};