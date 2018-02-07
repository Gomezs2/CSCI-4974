#include <iostream>
#include "product.h"

void Product::print() const {
	std::cout << "ID: " << id_ <<  "\tAverage-rating: " << avgRating_ << "\tASIN: " << asin_ << std::endl;
	std::cout << "Title: " << title_ << "\tGroup: " << group_  << std::endl;
	std::cout << "Similar Products: ";
	if (Simproducts_.size() == 0)
		std::cout << "None";
	for (int i = 0; i < Simproducts_.size(); ++i)
		std::cout << Simproducts_[i] << " ";
	std::cout << std::endl;
}