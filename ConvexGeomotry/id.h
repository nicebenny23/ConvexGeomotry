#pragma once
#include <stdint.h>
#include <limits>
namespace stn {
	static constexpr uint32_t none_id = std::numeric_limits<uint32_t>::max();



	template<typename T, typename Backing = std::uint32_t>
	struct typed_id {
		Backing id;
		constexpr explicit typed_id(Backing val) : id(val) {}
		constexpr bool operator==(const typed_id& other) const { return id == other.id; }
		constexpr bool operator!=(const typed_id& other) const { return id != other.id; }
		typed_id() = delete;
		constexpr explicit operator Backing() const { return id; }
		constexpr explicit operator size_t() const { return static_cast<size_t>(id); }
	};
	//nonmovable id type
	template<typename T>
	struct token_id {

		uint32_t id;
		constexpr explicit token_id(uint32_t val) : id(val) {}
		constexpr bool operator==(const token_id& other) const { return id == other.id; }
		constexpr bool operator!=(const token_id& other) const { return id != other.id; }
		token_id() = delete;
		token_id(token_id&& other) :id(other.id) {};
		token_id(const token_id& other) = delete;
		token_id& operator=(const token_id& other) = delete;
		token_id& operator=(token_id&& other) { id = other.id; return *this; };

		constexpr explicit operator uint32_t() const { return id; }
		constexpr explicit operator size_t() const { return static_cast<size_t>(id); }
	};

}
#include <functional>
template<typename T>
struct std::hash<stn::typed_id<T>> {
	size_t operator()(stn::typed_id<T> t) const noexcept {
		return std::hash<uint32_t>{}(t.id);
	}
};

template <typename T>
struct std::formatter<stn::typed_id<T>> : std::formatter<std::string> {
	template <typename FormatContext>
	auto format(const stn::typed_id<T>& v, FormatContext& ctx) const {
		return std::format_to(ctx.out(), "{}", v.id);
	}
};
