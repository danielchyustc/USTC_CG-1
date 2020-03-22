#pragma once

#include <Basic/TypeMap.h>

#include <Basic/Ptr.h>
#include <Engine/MeshEdit/AXAP.h>

#include <qobject.h>
#include <qtoolbox.h>

#include <map>

namespace Ubpa {
	class RawAPI_OGLW;
	class SObj;
	class Component;

	class Grid;

	enum WeightType {
		kEqual,
		kCotangent,
	};

	enum BoundShape {
		kCircle,
		kSquare,
	};

	enum ParaType {
		kASAP,
		kARAP,
	};

	class Attribute final {
	protected:
		Attribute();

	public:
		static Attribute* GetInstance() {
			static Attribute instance;
			return &instance;
		}

	public:
		void Init(QToolBox* tbox, RawAPI_OGLW* pOGLW);
		void SetSObj(Ptr<SObj> sobj);
		const Ptr<SObj> GetCurSObj() const { return curSObj.lock(); }
		void SetWeightType(WeightType weight_type) { weight_type_ = weight_type; }
		void SetBoundShape(BoundShape bound_shape) { bound_shape_ = bound_shape; }
		void SetParaType(ParaType para_type) { para_type_ = para_type; }
		template<typename T, typename = std::enable_if_t<std::is_base_of_v<Component, T>>>
		void SetCurCmpt() {
			auto target = componentType2item.find(typeid(T));
			if (target == componentType2item.end())
				return;

			tbox->setCurrentWidget(target->second);
		}

	private:
		void AddController(Ptr<SObj> sobj);

	private:
		class ComponentVisitor;
		friend class ComponentVisitor;

		QToolBox* tbox;

		TypeMap<QWidget*> componentType2item;
		std::map<QWidget*, Ptr<Grid>> item2grid;

		Ptr<ComponentVisitor> visitor;

		WPtr<SObj> curSObj;

		RawAPI_OGLW* pOGLW;

		WeightType weight_type_ = kEqual;
		BoundShape bound_shape_ = kCircle;
		ParaType para_type_ = kASAP;

		Ptr<AXAP> axap;
	};
}
