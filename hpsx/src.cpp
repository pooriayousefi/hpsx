
/***************************************************************************************/
/*                                                                                     */
/*                  high-performance descriptive statistics library                    */
/*                                                                                     */
/*                                    Pooria Yousefi                                   */
/*                                 pooriayousefi@aol.com                               */
/*                       https://www.linkedin.com/in/pooriayousefi/                    */
/*                                         2022                                        */
/*                                                                                     */
/***************************************************************************************/

// ------------------------------------------------
//
//               standard headers
//
// ------------------------------------------------

#include <type_traits>
#include <exception>
#include <execution>
#include <concepts>

#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>
#include <numeric>
#include <limits>
#include <random>
#include <cmath>

#include <initializer_list>
#include <valarray>
#include <vector>
#include <array>
#include <tuple>
#include <map>
#include <set>

// ------------------------------------------------
//
//             high-performance space
//
// ------------------------------------------------
namespace hp
{
	// --------------------------------------------
	//
	//                 constraints
	//
	// --------------------------------------------
	// arithmetic concept
	template<typename T>
	concept arithmetic = std::integral<T> || std::floating_point<T>;

	// integral value iterator concept
	template<typename T>
	concept integral_value_iterator = std::input_or_output_iterator<T> && std::integral<std::iter_value_t<T> >;

	// floating-point value iterator concept
	template<typename T>
	concept floating_point_value_iterator = std::input_or_output_iterator<T> && std::floating_point<std::iter_value_t<T> >;

	// arithmetic value iterator concept
	template<typename T>
	concept arithmetic_value_iterator = integral_value_iterator<T> || floating_point_value_iterator<T>;

	// arithmetic value iterator invocable concept = F f(ArgIter beg, ArgIter end)
	template<typename F, typename ArgIter>
	concept arithmetic_value_iterator_invocable = arithmetic_value_iterator<ArgIter> && std::invocable<F, ArgIter>;

	// --------------------------------------------
	//
	//                type definitions
	//
	// --------------------------------------------	
	// integrated numerical type
	template<typename T>
	using integrated_numerical_type = std::common_type_t<float, std::iter_value_t<T> >;

	// percentile typr
	template<arithmetic T>
	using percentile_t = std::array<T, 101Ui64>;

	// quartile type
	template<arithmetic T>
	using quartile_t = std::array<T, 5Ui64>;
}

// entry point
auto main()->int
{
	// --------------------------------------------
	//
	//               general utilities
	//
	// --------------------------------------------
	// summation
	auto Σ = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return std::reduce(
			std::execution::seq,
			beg,
			end,
			(T)0,
			std::plus<T>()
		);
	};

	// multiplication
	auto Π = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return std::reduce(
			std::execution::seq,
			beg,
			end,
			(T)1,
			std::multiplies<T>()
		);
	};

	// is even predicate
	auto is_even_value = []<std::integral T>(T value)
	{
		auto result{ false };
		if ((value % (T)2) == (T)0)
			result = true;
		return result;
	};

	// is odd predicate
	auto is_odd_value = []<std::integral T>(T value)
	{
		auto result{ false };
		if ((value % (T)2) != (T)0)
			result = true;
		return result;
	};

	// is all of the range values finite predicate
	auto is_all_of_the_values_finite = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		auto result{ false };
		if constexpr (std::is_integral_v<value_type>)
			result = true;
		if constexpr (std::is_floating_point_v<value_type>)
			if (std::all_of(beg, end, [](const auto& value) { return std::isfinite(value); }))
				result = true;
		return result;
	};

	// is all of the range values non-negative predicate
	auto is_all_of_the_values_nonnegative = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		auto result{ false };
		if (std::all_of(beg, end, [](const auto& value) { return value >= (value_type)0; }))
			result = true;
		return result;
	};

	// is all of the range values positive predicate
	auto is_all_of_the_values_positive = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		auto result{ false };
		if (std::all_of(beg, end, [](const auto& value) { return value > value_type(0); }))
			result = true;
		return result;
	};

	// is all of the range values non-positive predicate
	auto is_all_of_the_values_nonpositive = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		auto result{ false };
		if (std::all_of(beg, end, [](const auto& value) { return value <= value_type(0); }))
			result = true;
		return result;
	};

	// is all of the range values negative predicate
	auto is_all_of_the_values_negative = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		auto result{ false };
		if (std::all_of(beg, end, [](const auto& value) { return value < value_type(0); }))
			result = true;
		return result;
	};

	// histogram
	auto histogram = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using value_type = std::iter_value_t<Iter>;
		std::map<value_type, size_t, std::less<value_type> > hist{};
		auto it = beg;
		while (it != end)
		{
			++hist[*it];
			it = std::next(it);
		}
		return hist;
	};

	// contour
	auto contour = []<hp::arithmetic T>(const std::map<T, size_t, std::less<T> >& hist)
	{
		std::map<size_t, std::vector<T>, std::less<size_t> > invhist{};
		for (auto& elem : hist)
			invhist[elem.second].push_back(elem.first);
		return invhist;
	};

	// --------------------------------------------
	//
	//       1st-moments = central tendencies
	//
	// --------------------------------------------
	// Pythagorean arithmetic mean
	auto arithmetic_mean = [Σ]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		return Σ(beg, end) / std::distance(beg, end);
	};

	// Pythagorean geometric mean
	auto geometric_mean = [Π]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		return pow(Π(beg, end), 1.0f / std::distance(beg, end));
	};

	// Pythagorean harmonic mean
	auto harmonic_mean = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return std::distance(beg, end) /
			std::transform_reduce(
				std::execution::seq,
				beg,
				end,
				(T)0,
				std::plus<T>(),
				[](const auto& e) { return 1.0f / e; }
		);
	};

	// quadratic mean
	auto quadratic_mean = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return sqrt(std::transform_reduce(
			std::execution::seq,
			beg,
			end,
			(T)0,
			std::plus<T>(),
			[](const auto& e) { return e * e; }
		) / std::distance(beg, end));
	};

	// root mean square
	auto root_mean_square = quadratic_mean;

	// rms
	auto rms = root_mean_square;

	// cubic mean
	auto cubic_mean = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return cbrt(std::transform_reduce(
			std::execution::seq,
			beg,
			end,
			(T)0,
			std::plus<T>(),
			[](const auto& e) { return e * e * e; }
		) / std::distance(beg, end));
	};

	// generalized mean
	auto generalized_mean = []<std::floating_point P>(P p)
	{
		return[p]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
		{
			using T = hp::integrated_numerical_type<Iter>;
			return pow(std::transform_reduce(
				std::execution::seq,
				beg,
				end,
				(T)0,
				std::plus<T>(),
				[p](const auto& e) { return pow(e, p); }
			) / std::distance(beg, end), 1.0f / p);
		};
	};

	// power mean
	auto power_mean = generalized_mean;

	// Holder mean
	auto Holder_mean = power_mean;

	// contraharmonic mean
	auto contraharmonic_mean = [Σ]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		return std::transform_reduce(
			std::execution::seq,
			beg,
			end,
			(T)0,
			std::plus<T>(),
			[](const auto& e) { return e * e; }
		) / Σ(beg, end);
	};

	// Lehmer mean
	auto Lehmer_mean = []<std::floating_point P>(P p)
	{
		return[p]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
		{
			using T = hp::integrated_numerical_type<Iter>;
			return std::transform_reduce(
				std::execution::seq,
				beg,
				end,
				(T)0,
				std::plus<T>(),
				[p](const auto& e) { return pow(e, p); }
			) / std::transform_reduce(
				std::execution::seq,
				beg,
				end,
				(T)0,
				std::plus<T>(),
				[p](const auto& e) {return pow(e, p - 1.0f); }
			);
		};
	};

	// weighted Pythagorean arithmetic mean
	auto weighted_arithmetic_mean = [Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		auto Σw{ Σ(wtbeg, wtend) };
		return[Σw, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * x; }
			) / Σw;
		};
	};

	// weighted Pythagorean geometric mean
	auto weighted_geometric_mean = [Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		auto Σw{ Σ(wtbeg, wtend) };
		return[Σw, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return exp(std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * log((T)x); }
			) / Σw);
		};
	};

	// weighted Pythagorean harmonic mean
	auto weighted_harmonic_mean = [Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		auto Σw{ Σ(wtbeg, wtend) };
		return[Σw, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return Σw / std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w / (T)x; });
		};
	};

	// weighted quadratic mean
	auto weighted_quadratic_mean = [Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		auto Σw{ Σ(wtbeg, wtend) };
		return[Σw, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return sqrt(std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * x * x; }
			) / Σw);
		};
	};

	// weighted root mean square
	auto weighted_root_mean_square = weighted_quadratic_mean;

	// weighted rms
	auto weighted_rms = weighted_root_mean_square;

	// weighted cubic mean
	auto weighted_cubic_mean = [Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		auto Σw{ Σ(wtbeg, wtend) };
		return[Σw, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return cbrt(std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * x * x * x; }
			) / Σw);
		};
	};

	// weighted generalized mean
	auto weighted_generalized_mean = [Σ]<std::floating_point P>(P p)
	{
		return[p, Σ]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
		{
			using W = hp::integrated_numerical_type<WtIter>;
			auto Σw{ Σ(wtbeg, wtend) };
			return[Σw, wtbeg, wtend, p]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
			{
				using X = hp::integrated_numerical_type<XIter>;
				using T = std::common_type_t<W, X>;
				return pow(std::transform_reduce(
					std::execution::seq,
					wtbeg,
					wtend,
					xbeg,
					(T)0,
					std::plus<T>(),
					[p](const auto& w, const auto& x) { return w * pow(x, p); }
				) / Σw, 1.0f / p);
			};
		};
	};

	// weighted power mean
	auto weighted_power_mean = weighted_generalized_mean;

	// weighted Holder mean
	auto weighted_Holder_mean = weighted_power_mean;

	// weighted contraharmonic mean
	auto weighted_contraharmonic_mean = []<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
	{
		using W = hp::integrated_numerical_type<WtIter>;
		return[wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			using X = hp::integrated_numerical_type<XIter>;
			using T = std::common_type_t<W, X>;
			return std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * x * x; }
			) / std::transform_reduce(
				std::execution::seq,
				wtbeg,
				wtend,
				xbeg,
				(T)0,
				std::plus<T>(),
				[](const auto& w, const auto& x) { return w * x; }
			);
		};
	};

	// weighted Lehmer mean
	auto weighted_Lehmer_mean = []<std::floating_point P>(P p)
	{
		return [p]<hp::arithmetic_value_iterator WtIter>(WtIter wtbeg, WtIter wtend)
		{
			using W = hp::integrated_numerical_type<WtIter>;
			return[p, wtbeg, wtend]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
			{
				using X = hp::integrated_numerical_type<XIter>;
				using T = std::common_type_t<W, X>;
				return std::transform_reduce(
					std::execution::seq,
					wtbeg,
					wtend,
					xbeg,
					(T)0,
					std::plus<T>(),
					[p](const auto& w, const auto& x) { return w * pow(x, p); }
				) / std::transform_reduce(
					std::execution::seq,
					wtbeg,
					wtend,
					xbeg,
					(T)0,
					std::plus<T>(),
					[p](const auto& w, const auto& x) { return w * pow(x, p - 1.0f); }
				);
			};
		};
	};

	// kth-moment
	auto kth_moment = []<std::floating_point K>(K k)
	{
		return [k]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			return[k, xbeg, xend]<hp::arithmetic_value_iterator XKIter>(XKIter xkbeg, XKIter xkend)
			{
				std::transform(
					std::execution::par_unseq,
					xbeg,
					xend,
					xkbeg,
					[k](const auto& x) { return pow(x, k); });
				return[xkbeg, xkend]<hp::arithmetic_value_iterator_invocable F>(F EV)
				{
					return std::invoke(EV, xkbeg, xkend);
				};
			};
		};
	};

	// kth central moment
	auto kth_central_moment = []<std::floating_point K>(K k)
	{
		return [k]<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
		{
			return[k, xbeg, xend]<hp::arithmetic_value_iterator XKIter>(XKIter xkbeg, XKIter xkend)
			{
				return[k, xbeg, xend, xkbeg, xkend]<hp::arithmetic_value_iterator_invocable F>(F EV)
				{
					auto μ{ std::invoke(EV, xbeg, xend) };
					std::transform(
						std::execution::par_unseq,
						xbeg,
						xend,
						xkbeg,
						[μ, k](const auto& x) { return pow(x - μ, k); }
					);
					return std::invoke(EV, xkbeg, xkend);
				};
			};
		};
	};

	// median
	// NOTE! range must be sorted otherwise, median value is not correct.
	auto median = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		auto n{ std::distance(beg, end) };
		auto result{ (T)0 };
		if ((n % 2Ui64) != 0Ui64) {
			auto it = beg;
			for (auto i{ 1Ui64 }; i <= n / 2Ui64; ++i)
				it = std::next(it);
			result = (T)*it;
		}
		else {
			auto it = beg;
			for (auto i{ 1Ui64 }; i <= n / 2Ui64 - 1Ui64; ++i)
				it = std::next(it);
			auto v{ *it };
			it = std::next(it);
			result = (T)(v + *it) / 2.0f;
		}
		return result;
	};

	// minimum/maximum/range/mid-range
	auto min_max_range_midrange = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		auto [minval_iter, maxval_iter] = std::minmax_element(beg, end);
		return std::make_tuple(
			*minval_iter,
			*maxval_iter,
			*maxval_iter - *minval_iter,
			(*maxval_iter + *minval_iter) / 2.0f
		);
	};

	// minimum/maximum/range/mid-extreme
	auto min_max_range_midextreme = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		auto [minval_iter, maxval_iter] = std::minmax_element(beg, end);
		return std::make_tuple(
			*minval_iter,
			*maxval_iter,
			*maxval_iter - *minval_iter,
			(*maxval_iter + *minval_iter) / 2.0f
		);
	};

	// mode
	auto mode = [histogram, contour]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		auto hist{ histogram(beg, end) };
		auto cont{ contour(hist) };
		return std::prev(cont.end())->second;
	};

	// ------------------------------------------------
	//
	//                  dispersions
	//
	// ------------------------------------------------
	// mean absolute deviation (MAD) around the mean (μ)
	auto mean_absolute_deviation_around_the_mean = []<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		auto n{ std::distance(beg, end) };
		return[beg, end, n]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			auto μ{ std::invoke(EV, beg, end) };
			return std::transform_reduce(
				std::execution::seq,
				beg,
				end,
				(T)0,
				std::plus<T>(),
				[μ](const auto& x) { return abs(x - μ); }
			) / n;
		};
	};
	auto average_absolute_deviation_around_the_mean = mean_absolute_deviation_around_the_mean;

	// mean absolute deviation (MAD) around the median
	auto mean_absolute_deviation_around_the_median = [median]<hp::arithmetic_value_iterator Iter>(Iter beg, Iter end)
	{
		using T = hp::integrated_numerical_type<Iter>;
		auto med{ median(beg, end) };
		return std::transform_reduce(
			std::execution::seq,
			beg,
			end,
			(T)0,
			std::plus<T>(),
			[med](const auto& x) { return abs(x - med); }
		) / std::distance(beg, end);
	};
	auto average_absolute_deviation_around_the_median = mean_absolute_deviation_around_the_median;

	// variance	
	auto variance = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			return std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ);
		};
	};

	// standard deviation
	auto standard_deviation = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			return sqrt(std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ));
		};
	};

	// coefficient of variation
	auto coefficient_of_variation = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto σ{ sqrt(std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ)) };
			return σ / μ;
		};
	};
	auto CV = coefficient_of_variation;
	auto CoV = CV;
	auto relative_standard_deviation = CoV;
	auto RSD = relative_standard_deviation;

	// index of dispersion
	auto index_of_dispersion = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto σ2{ std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ) };
			return σ2 / μ;
		};
	};
	auto dispersion_index = index_of_dispersion;
	auto coefficient_of_dispersion = dispersion_index;
	auto relative_variance = coefficient_of_dispersion;
	auto variance_to_mean_ratio = relative_variance;
	auto VMR = variance_to_mean_ratio;
	auto Fano_factor = VMR;

	// signal to noise ratio (SNR)
	auto signal_to_noise_ratio = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto σ{ sqrt(std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ)) };
			return μ / σ;
		};
	};
	auto SNR = signal_to_noise_ratio;

	// efficiency
	auto efficiency = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto σ2{ std::invoke(EV, std::begin(x2), std::end(x2)) - (μ * μ) };
			return σ2 / (μ * μ);
		};
	};

	// skewness
	auto skew = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto m2{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto σ2{ m2 - (μ * μ) };
			std::for_each(
				std::execution::par_unseq,
				std::begin(x2),
				std::end(x2),
				[μ](auto& x)
				{
					x = sqrt(x);
					x = pow(x - μ, 3.0f);
				}
			);
			auto μ3{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			return μ3 / pow(σ2, 1.5f);
		};
	};

	// kurtosis
	auto kurtosis = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto m2{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto σ2{ m2 - (μ * μ) };
			std::for_each(
				std::execution::par_unseq,
				std::begin(x2),
				std::end(x2),
				[μ](auto& x)
				{
					x = sqrt(x);
					x = pow(x - μ, 4.0f);
				}
			);
			auto μ4{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			return μ4 / (σ2 * σ2);
		};
	};

	// excess kurtosis
	auto excess_kurtosis = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };
			auto m2{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto σ2{ m2 - (μ * μ) };
			std::for_each(
				std::execution::par_unseq,
				std::begin(x2),
				std::end(x2),
				[μ](auto& x)
				{
					x = sqrt(x);
					x = pow(x - μ, 4.0f);
				}
			);
			auto μ4{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			return μ4 / (σ2 * σ2) - 3.0f;
		};
	};

	// dispersion
	auto dispersion = []<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return[xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
		{
			using T = hp::integrated_numerical_type<Iter>;
			std::vector<T> x2(std::distance(xbeg, xend), (T)0);
			std::transform(
				std::execution::seq,
				xbeg,
				xend,
				std::begin(x2),
				[](const auto& x) { return x * x; }
			);
			auto μ{ std::invoke(EV, xbeg, xend) };// ------------------------------------ mean
			auto m2{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto σ2{ m2 - (μ * μ) };// -------------------------------------------------- variance
			auto σ{ sqrt(σ2) };// ------------------------------------------------------- standard deviation
			auto cov{ σ / μ };// -------------------------------------------------------- coefficient of variation
			auto vmr{ σ2 / μ };// ------------------------------------------------------- coefficient of dispersion (VMR)
			auto snr{ μ / σ };// -------------------------------------------------------- signal to noise ratio (SNR)
			auto eff{ σ2 / (μ * μ) };// ------------------------------------------------- efficiency
			std::for_each(
				std::execution::par_unseq,
				std::begin(x2),
				std::end(x2),
				[μ](auto& x)
				{
					x = sqrt(x);
					x = pow(x - μ, 3.0f);
				}
			);
			auto μ3{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto skw{ μ3 / pow(σ2, 1.5f) };// ------------------------------------------- skewness
			std::for_each(
				std::execution::par_unseq,
				std::begin(x2),
				std::end(x2),
				[μ](auto& x)
				{
					x = cbrt(x) + μ;
					x = pow(x - μ, 4.0f);
				}
			);
			auto μ4{ std::invoke(EV, std::begin(x2), std::end(x2)) };
			auto xkurt{ μ4 / (σ2 * σ2) - 3.0f };// -------------------------------------- excess kurtosis
			return std::make_tuple(μ, σ2, σ, cov, vmr, snr, eff, skw, xkurt);
		};
	};
	auto variability = dispersion;

	// covariance
	auto covariance = []<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
	{
		using X = hp::integrated_numerical_type<XIter>;
		return[xbeg, xend]<hp::arithmetic_value_iterator YIter>(YIter ybeg, YIter yend)
		{
			using Y = hp::integrated_numerical_type<YIter>;
			return[xbeg, xend, ybeg, yend]<hp::arithmetic_value_iterator_invocable F>(F EV)
			{
				using T = std::common_type_t<X, Y>;
				auto μx{ std::invoke(EV, xbeg, xend) };
				auto μy{ std::invoke(EV, ybeg, yend) };
				
				std::vector<T> xy(std::distance(xbeg, xend), (T)0);
				std::transform(
					std::execution::seq,
					xbeg,
					xend,
					ybeg,
					std::begin(xy),
					[](const auto& xi, const auto& yi) { return xi * yi; }
				);
				return std::invoke(EV, std::begin(xy), std::end(xy)) - (μx * μy);
			};
		};
	};
	auto cov = covariance;
	auto COV = cov;

	// Pearson correlation coefficient (PCC)
	auto correlation_coefficient = []<hp::arithmetic_value_iterator XIter>(XIter xbeg, XIter xend)
	{
		using X = hp::integrated_numerical_type<XIter>;
		return[xbeg, xend]<hp::arithmetic_value_iterator YIter>(YIter ybeg, YIter yend)
		{
			using Y = hp::integrated_numerical_type<YIter>;
			return[xbeg, xend, ybeg, yend]<hp::arithmetic_value_iterator_invocable F>(F EV)
			{
				using T = std::common_type_t<X, Y>;
				auto μx{ std::invoke(EV, xbeg, xend) };
				auto μy{ std::invoke(EV, ybeg, yend) };
				std::vector<T> x2(std::distance(xbeg, xend), (T)0);
				std::transform(
					std::execution::seq,
					xbeg,
					xend,
					std::begin(x2),
					[](const auto& x) { return x * x; }
				);
				std::vector<T> y2(std::distance(ybeg, yend), (T)0);
				std::transform(
					std::execution::seq,
					ybeg,
					yend,
					std::begin(y2),
					[](const auto& y) { return y * y; }
				);
				std::vector<T> xy(std::distance(xbeg, xend), (T)0);
				std::transform(
					std::execution::seq,
					xbeg,
					xend,
					ybeg,
					std::begin(xy),
					[](const auto& xi, const auto& yi) { return xi * yi; }
				);
				auto σ2x{ std::invoke(EV, std::begin(x2), std::end(x2)) - (μx * μx) };
				auto σ2y{ std::invoke(EV, std::begin(y2), std::end(y2)) - (μy * μy) };
				auto covxy{ std::invoke(EV, std::begin(xy), std::end(xy)) - (μx * μy) };				
				return covxy / sqrt(σ2x * σ2y);
			};
		};
	};
	auto Pearson_correlation_coefficient = correlation_coefficient;
	auto PCC = Pearson_correlation_coefficient;
	auto Pearsons_r = PCC;
	auto Pearson_product_moment_correlation_coefficient = Pearsons_r;
	auto PPMCC = Pearson_product_moment_correlation_coefficient;
	auto bivariate_correlation = PPMCC;
	auto ρxy = bivariate_correlation;

	// moment generating function (MGF)
	auto moment_generating_function = []<std::floating_point T>(T t)
	{
		return[t]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
		{
			using X = hp::integrated_numerical_type<Iter>;
			return[t, xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
			{
				using Y = std::common_type_t<T, X>;
				std::vector<Y> y(std::distance(xbeg, xend), (Y)0);
				std::transform(
					std::execution::seq,
					xbeg, 
					xend,
					std::begin(y),
					[t](const auto& x) { return exp(t * x); }
				);
				return std::invoke(EV, std::begin(y), std::end(y));
			};
		};
	};
	auto MGF = moment_generating_function;

	// cummulant generating function
	auto cummulant_generating_function = []<std::floating_point T>(T t)
	{
		return[t]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
		{
			using X = hp::integrated_numerical_type<Iter>;
			return[t, xbeg, xend]<hp::arithmetic_value_iterator_invocable F>(F EV)
			{
				using Y = std::common_type_t<T, X>;
				std::vector<Y> y(std::distance(xbeg, xend), (Y)0);
				std::transform(
					std::execution::seq,
					xbeg,
					xend,
					std::begin(y),
					[t](const auto& x) { return exp(t * x); }
				);
				return log(std::invoke(EV, std::begin(y), std::end(y)));
			};
		};
	};
	auto K = cummulant_generating_function;

	// ------------------------------------------------
	//
	//                   percentile
	//
	// ------------------------------------------------
	// input range in this function must be incrementally sorted
	auto percentile = [is_odd_value]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		using T = hp::integrated_numerical_type<Iter>;
		hp::percentile_t<T> p{};
		size_t N{ std::distance(xbeg, xend) };
		p[0Ui64] = *xbeg;
		p[100Ui64] = *std::prev(xend);
		if (is_odd_value(N))
			for (auto i{ 1Ui64 }; i < 100Ui64; ++i)
				p[i] = *std::next(xbeg, i * N / 100Ui64);
		else
			for (auto i{ 1Ui64 }; i < 100Ui64; ++i)
			{
				auto v1{ *std::next(xbeg, i * N / 100Ui64 - 1Ui64) };
				auto v2{ *std::next(xbeg, i * N / 100Ui64) };
				p[i] = ((v1 + v2) / (T)2);
			}
		return p;
	};

	// ------------------------------------------------
	//
	//                    quartile
	//
	// ------------------------------------------------
	// input range in this function must be incrementally sorted
	auto quartile = [is_odd_value]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		using T = hp::integrated_numerical_type<Iter>;
		hp::quartile_t<T> Q{};
		auto N{ std::distance(xbeg, xend) };
		Q[0] = *xbeg;
		Q[4] = *std::prev(xend);
		if (is_odd_value(N))
		{
			Q[1] = *std::next(xbeg, N / 4Ui64);
			Q[2] = *std::next(xbeg, N / 2Ui64);
			Q[3] = *std::next(xbeg, 3Ui64 * N / 4Ui64);
		}
		else
		{
			Q[1] = (*std::next(xbeg, N / 4Ui64 - 1Ui64) + *std::next(xbeg, N / 4Ui64)) / (T)2;
			Q[2] = (*std::next(xbeg, N / 2Ui64 - 1Ui64) + *std::next(xbeg, N / 2Ui64)) / (T)2;
			Q[3] = (*std::next(xbeg, 3Ui64 * N / 4Ui64 - 1Ui64) + *std::next(xbeg, 3Ui64 * N / 4Ui64)) / (T)2;
		}
		return Q;
	};

	// inter-quartile range
	auto inter_quartile_range = [quartile]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		auto Q = quartile(xbeg, xend);
		return Q[3] - Q[1];
	};
	auto mid_spread = inter_quartile_range;

	// quartile deviation
	auto quartile_deviation = [inter_quartile_range]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		return inter_quartile_range(xbeg, xend) / 2.0f;
	};

	// mid-hinge
	auto mid_hinge = [quartile]<hp::arithmetic_value_iterator Iter>(Iter xbeg, Iter xend)
	{
		auto Q = quartile(xbeg, xend);
		return (Q[3] + Q[1]) / 2.0f;
	};

	try
	{
		// test
		auto test = [&]()
		{
			// general utilities test
			auto general_utilities_test = [&]()
			{
				std::vector<float> rng0{ 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
				auto b{ std::begin(rng0) };
				auto e{ std::end(rng0) };
				std::cout
					<< "\nrange: " << [&rng0]()
				{
					for (const auto& elem : rng0)
						std::cout << elem << ' ';
					return ' ';
				}()
					<< "\nrange sum: " << Σ(b, e)
					<< "\nrange multiplication: " << Π(b, e)
					<< "\nis all of the values finite: " << [b, e, is_all_of_the_values_finite]()
				{
					std::string s{ "" };
					if (is_all_of_the_values_finite(b, e))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis all of the values nonnegative: " << [b, e, is_all_of_the_values_nonnegative]()
				{
					std::string s{ "" };
					if (is_all_of_the_values_nonnegative(b, e))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis all of the values positive: " << [b, e, is_all_of_the_values_positive]()
				{
					std::string s{ "" };
					if (is_all_of_the_values_positive(b, e))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis all of the values nonpositive: " << [b, e, is_all_of_the_values_nonpositive]()
				{
					std::string s{ "" };
					if (is_all_of_the_values_nonpositive(b, e))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis all of the values negative: " << [b, e, is_all_of_the_values_negative]()
				{
					std::string s{ "" };
					if (is_all_of_the_values_negative(b, e))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis -2 odd: " << [is_odd_value]()
				{
					std::string s{ "" };
					if (is_odd_value(-2))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis 3U odd: " << [is_odd_value]()
				{
					std::string s{ "" };
					if (is_odd_value(3U))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis 3Ui64 even: " << [is_even_value]()
				{
					std::string s{ "" };
					if (is_even_value(3Ui64))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< "\nis -88i64 even: " << [is_even_value]()
				{
					std::string s{ "" };
					if (is_even_value(-88i64))
						s += "yes";
					else
						s += "no";
					return s;
				}()
					<< std::endl;
				std::multiset<int> rng1{ -2, -2, 0, 3, -2, 1, 0, -2, 3, 4, 0 };
				std::cout
					<< "\n\nanother range: " << [&rng1]()
				{
					for (const auto& elem : rng1)
						std::cout << elem << ' ';
					return ' ';
				}()
					<< "\nhistogram: " << [&rng1, histogram, contour]()
				{
					auto b{ std::begin(rng1) };
					auto e{ std::end(rng1) };
					auto hist{ histogram(b, e) };
					for (const auto& elem : hist)
						std::cout << '(' << elem.first << ',' << elem.second << ") ";
					auto cont{ contour(hist) };
					std::cout << "\ncontour: ";
					for (const auto& elem : cont)
						std::cout << '(' << elem.first << ',' << [&elem]()
					{
						for (const auto& elems : elem.second)
							std::cout << elems << ' ';
						return ' ';
					}() << ") ";
					return ' ';
				}()
					<< std::endl;
			};
			general_utilities_test();

			// central tendencies test
			auto central_tendencies_test = [&]()
			{

			};
			central_tendencies_test();

			// dispersions test
			auto dispersions_test = [&]()
			{

			};
			dispersions_test();

			// percentile test
			auto percentile_test = [&]()
			{

			};
			percentile_test();

			// quartile test
			auto quartile_test = [&]()
			{

			};
			quartile_test();
		};
		test();		
		
		return EXIT_SUCCESS;
	}
	catch (const std::exception& xxx)
	{
		std::cout << xxx.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cout << "UNCAUGHT EXCEPTION DETECTED" << std::endl;
		return EXIT_FAILURE;
	}
}