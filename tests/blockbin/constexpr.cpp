//  constexpr.cpp : compile-time tests for constexpr of blockbinary type
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <strstream>

//#include <universal/integer/integer.hpp>
#include <universal/blockbin/blockbinary.hpp>

int main()
try {
	using namespace std;
	using namespace sw::unum;

	std::string tag = "blockbinary storage class constexpr compile-time testing";

	{
		constexpr blockbinary<8> b8;
		constexpr blockbinary<8, uint8_t> b8_1b(0x5555);
		constexpr blockbinary<8, uint16_t> b8_2b(0x5555);
		constexpr blockbinary<8, uint32_t> b8_4b(0x5555);

		stringstream ss;
		ss << b8 << '\n' << b8_1b << '\n' << b8_2b << '\n' << b8_4b << endl;
		cout << ss.str() << endl;
	}

	{
		constexpr blockbinary<16> b16;
		constexpr blockbinary<16, uint8_t> b16_2(0x5555);
		constexpr blockbinary<16, uint16_t> b16_1(0x5555);
		constexpr blockbinary<16, uint32_t> b16_4b(0x5555);

		stringstream ss;
		ss << b16 << '\n' << b16_1 << '\n' << b16_2 << '\n' << b16_4b << endl;
		cout << ss.str() << endl;
	}

}
catch (char const* msg) {
	std::cerr << msg << '\n';
	return EXIT_FAILURE;
}
catch (const std::runtime_error& err) {
	std::cerr << "Uncaught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << '\n';
	return EXIT_FAILURE;
}